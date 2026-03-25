#' Inverse Distance Weight in graphs
#'
#'@description
#'Oblicza idw gdzie waga to odleglosc mierzona wzdluz drog
#'a nie nakrotsza droga jak w klasycznym idw.
#'
#'@details
#'Waga to dlugosc odcinka drogi razy liczba odcinkow
#'
#' @param roads Obiekt klasy sf z siecia drog.
#' @param aprs Data.frame o strukturze: 2 wiersze, 3 kolumny
#'        W kolumnach apry
#'        W pierwszym wierszu współrzędne apru w formacie "20.025 52.60337"
#'        W drugim wierszu wartość r dla apru w formacie "58.11"
#'        Wszystko jako tekst.
#' @param power
#'
#' @return sf object
#' @export
#'
#' @examples
#' test <- idw_graph(roads, aprs)
#'
idw_graph <- function(roads, aprs, power = 2) {
  library(sf)
  library(dplyr)

  # Trzeba usunac duplikaty drog, ktore sie roznia tylko kierunkiem
  tictoc::tic()
  roads <- roads %>%
    mutate(centr = st_centroid(geometry)) %>%
    distinct(centr, .keep_all = TRUE) %>%
    select(-centr)
  aprs1 <- stringr::str_trim(aprs[1, ])
  tictoc::toc()

  ## APRy do mapy
  to_points <- function(x) {
    lon <- as.numeric(stringr::str_split(x[1], " ")[[1]][1])
    lat <- as.numeric(stringr::str_split(x[1], " ")[[1]][2])
    data.frame(lon = lon, lat = lat)
  }

  pprs <- purrr::map(aprs1, to_points) %>% data.table::rbindlist()
  pprs <- st_as_sf(pprs, coords = c("lon", "lat"), crs = st_crs(roads))
  pprs$r <- as.numeric(t(aprs[2, ]))

  ## test wizualny
  # library(tmap)
  # tmap_mode("view")
  #
  # tm_shape(roads) + tm_lines() + tm_shape(pprs) + tm_dots(col = "red")
  # punkty nie pasuja do drog

  ## dopasowanie pprow do drog i przypisanie r
  roads$r <- NA

  tictoc::tic()
  for (i in 1:nrow(pprs)) {
    points_align <- st_nearest_feature(pprs[i, ], roads)
    roads$r[points_align] <- pprs$r[i]
  }
  tictoc::toc()

  #weryfikacja dopasowania
  # tm_shape(filter(roads, !is.na(r))) + tm_lines() + tm_shape(pprs) + tm_dots(col = "red")

  # Tworzenie nodes
  n <- roads %>% select(r, dl_segm) %>% mutate(nodeID = c(1:n()))

  # Tworzenie edges miedzy nodes
  tictoc::tic()
  e <- n %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(nodeID = L1) %>%
    group_by(nodeID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('a', 'b'), times = n() / 2))
  tictoc::toc()
  # Otrzymujemy dwa edges na jeden node - a i b
  # trzeba je pogrupowac wg xy wtedy jeden bedzie from a drugi to. Te bez pary trzeba wyrzucic

  tictoc::tic()
  e <- e %>% mutate(xy = paste(.$X, .$Y)) %>% group_by(xy)
  e$EdgeID <- group_indices(e)
  e <- e %>% group_by(EdgeID) %>% mutate(n = n()) %>% filter(n > 1) # n = 1 oznacza slepy koniec drogi wiec  sie nie nadaje na edge
  tictoc::toc()

  # Jest problem ze skrzyzowaniami
  #tm_shape(filter(n, nodeID %in% c(1, 19, 247, 327))) + tm_lines(col = "red")

  # Trzeba potworzyc nowe edge na zasadzie kazdy node z kazdym (expand.grid)
  tictoc::tic()
  e2 <- filter(e, n > 2)
  e_list <- collapse::rsplit(e2, e2$EdgeID)
  krzyz <- function(x) {
    expand.grid(x$nodeID, x$nodeID) %>%
      filter(Var1 != Var2) %>%
      mutate(EdgeID = unique(x$EdgeID))
  }
  e2 <- e_list %>%
    purrr::map(\(x) krzyz(x)) %>%
    data.table::rbindlist() %>%
    rename(all_of(c(from = "Var1", to = "Var2")))
  e1 <- filter(e, n == 2) %>%
    select(nodeID, EdgeID) %>%
    mutate(start_end = rep(c('from', 'to'), times = n() / 2)) %>%
    tidyr::pivot_wider(values_from = nodeID, names_from = start_end)

  e <- rbind(e1, e2)
  tictoc::toc()

  tictoc::tic()

  # edges to linie wiec trzeba dodac dlugosc, ktora bedzie waga dla shortest_paths
  # ponieważ n1 to centroida odcinka to długość to suma połowy odcinków z obu stron, czyli 0.5*dl_odcinka1 + 0.5*dl_odcinka2.
  # Dolacz geometrie wezlow "from" i "to" do krawedzi i policz odleglosci wektorowo
  n <- roads %>% select(r, dl_segm) %>% mutate(nodeID = c(1:n()))
  n1_lookup <- n %>% st_drop_geometry() |> dplyr::select(nodeID, dl_segm)
  e_dist <- e %>%
    dplyr::left_join(n1_lookup, by = c("from" = "nodeID")) %>%
    dplyr::rename(dl_from = dl_segm / 2) %>%
    dplyr::left_join(n1_lookup, by = c("to" = "nodeID")) %>%
    dplyr::rename(dl_to = dl_segm / 2) %>%
    dplyr::mutate(
      waga = dl_from + dl_to # Można dodać jeszcze poprawkę na prędkość.
    ) |>
    select(EdgeID, from, to, waga)
  tictoc::toc()

  library(tidygraph)
  # Tworzenie i poprawianie grafu
  graph <- tbl_graph(
    nodes = as_tibble(n),
    edges = as_tibble(e_dist),
    directed = FALSE
  )

  # IDW z pprow wg drog ----
  ppr_node <- graph |>
    activate(nodes) |>
    as_tibble() |>
    dplyr::filter(!is.na(r))
  src_nodes <- ppr_node$nodeID
  src_vals <- ppr_node$r
  K <- length(src_nodes)
  odl_names <- paste0("odl", seq_len(K))

  # Obliczenia odleglosci - macierz Dijkstra bez sf w workersach
  tictoc::tic("Dijkstra distances")

  ig <- igraph::graph_from_data_frame(
    d = as_tibble(e_dist) |> dplyr::select(from, to, waga),
    directed = FALSE,
    vertices = as_tibble(n) |> sf::st_drop_geometry() |> dplyr::select(nodeID)
  )

  dist_mat <- igraph::distances(
    ig,
    v = as.character(src_nodes),
    to = igraph::V(ig),
    mode = "all",
    weights = igraph::E(ig)$waga,
    algorithm = "dijkstra"
  )
  # dist_mat jest K x N, potrzebujemy N x K
  dist_mat <- t(dist_mat)
  colnames(dist_mat) <- odl_names

  tictoc::toc()

  # Obliczanie IDW (wektorowo, bez petli)
  tictoc::tic("IDW")
  dist_mat <- as.matrix(dist_mat)

  zero_mask <- dist_mat == 0
  has_zero <- apply(zero_mask, 1, any, na.rm = TRUE)

  w <- 1 / (dist_mat^power)
  w[!is.finite(w)] <- 0

  idw_est <- as.vector((w %*% src_vals) / rowSums(w))

  if (any(has_zero)) {
    src_idx <- max.col(zero_mask * 1L, ties.method = "first")
    idw_est[has_zero] <- src_vals[src_idx[has_zero]]
  }

  nodes_tbl <- as_tibble(n) |> sf::st_drop_geometry()
  wynik <- n |>
    dplyr::mutate(r = dplyr::if_else(is.na(r), idw_est, r))
  tictoc::toc()

  return(wynik)
}
