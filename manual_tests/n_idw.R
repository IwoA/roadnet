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
#' test <- n_idw(roads, aprs)
#'
n_idw <- function(roads, aprs, power = 2) {
  library(sf)
  library(dplyr)

  # Trzeba usunac duplikaty drog, ktore sie roznia tylko kierunkiem
  roads <- roads %>%
    mutate(centr = st_centroid(geometry)) %>%
    distinct(centr, .keep_all = TRUE) %>%
    select(-centr)
  aprs1 <- stringr::str_trim(aprs[1, ])

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

  for (i in 1:nrow(pprs)) {
    points_align <- st_nearest_feature(pprs[i, ], roads)
    roads$r[points_align] <- pprs$r[i]
  }

  #weryfikacja dopasowania
  # tm_shape(filter(roads, !is.na(r))) + tm_lines() + tm_shape(pprs) + tm_dots(col = "red")

  # Tworzenie nodes
  n <- roads %>% select(r) %>% mutate(nodeID = c(1:n()))

  # Tworzenie edges miedzy nodes
  e <- n %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(nodeID = L1) %>%
    group_by(nodeID) %>%
    slice(c(1, n())) %>%
    ungroup() %>%
    mutate(start_end = rep(c('a', 'b'), times = n() / 2))
  # Otrzymujemy dwa edges na jeden node - a i b
  # trzeba je pogrupowac wg xy wtedy jeden bedzie from a drugi to. Te bez pary trzeba wyrzucic

  e <- e %>% mutate(xy = paste(.$X, .$Y)) %>% group_by(xy)
  e$EdgeID <- group_indices(e)
  e <- e %>% group_by(EdgeID) %>% mutate(n = n()) %>% filter(n > 1) # n = 1 oznacza slepy koniec drogi wiec  sie nie nadaje na edge

  # Jest problem ze skrzyzowaniami
  #tm_shape(filter(n, nodeID %in% c(1, 19, 247, 327))) + tm_lines(col = "red")

  # Trzeba potworzyc nowe edge na zasadzie kazdy node z kazdym (expand.grid)
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

  n1 <- sf::st_centroid(n) # nodes to punkty. To musi byc tu, zeby wagi e sie policzyly
  # edges to linie wiec trzeba dodac dlugosc, ktora bedzie waga dla shortest_paths
  # tictoc::tic()
  e_dist <- e %>%
    rowwise() %>%
    mutate(
      waga = round(
        sf::st_distance(
          filter(n1, nodeID == from),
          filter(n1, nodeID == to)
        ),
        3
      )
    )
  # tictoc::toc()

  library(tidygraph)
  # Tworzenie i poprawianie grafu
  graph <- tbl_graph(
    nodes = as_tibble(n),
    edges = as_tibble(e_dist),
    directed = FALSE
  )

  # IDW z pprow wg drog ----
  ppr_node <- graph %>% activate(nodes) %>% as_tibble() %>% filter(!is.na(r))
  new_n <- select(n, nodeID, r, geometry) %>%
    mutate(odl1 = NA, odl2 = NA, odl3 = NA)
  # Na razie jako nowy graf
  graph1 <- tbl_graph(
    nodes = new_n,
    edges = as_tibble(e_dist),
    directed = FALSE
  )

  #tictoc::tic()
  for (i in 1:nrow(new_n)) {
    # Obliczenia odleglosci kazdej drogi od kazdego z 3 pprow
    cat("wiersz ", i, "z ", nrow(new_n), "\r")

    graph1 <- graph1 %>%
      mutate(
        odl1 = ifelse(
          nodeID == new_n$nodeID[i],
          node_distance_from(
            ppr_node$nodeID[1],
            algorithm = "dijkstra",
            weights = waga
          ),
          odl1
        ),
        odl2 = ifelse(
          nodeID == new_n$nodeID[i],
          node_distance_from(
            ppr_node$nodeID[2],
            algorithm = "dijkstra",
            weights = waga
          ),
          odl2
        ),
        odl3 = ifelse(
          nodeID == new_n$nodeID[i],
          node_distance_from(
            ppr_node$nodeID[3],
            algorithm = "dijkstra",
            weights = waga
          ),
          odl3
        ),
      )
  }
  #tictoc::toc()

  # Obliczanie idw
  wynik <- graph1 %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(
      odl1 = 1 / odl1^power,
      odl2 = 1 / odl2^power,
      odl3 = 1 / odl3^power
    ) %>% # wagi - inverse distance
    tidyr::pivot_longer(starts_with("odl"), names_to = "odl") %>% # zmiana formatu do obliczenia calkowitej odl
    group_by(nodeID) %>%
    mutate(suma_odl = sum(value)) %>% # suma odleglosci
    tidyr::pivot_wider(names_from = odl, values_from = value) %>% # powrot formatu
    mutate(
      odl1 = odl1 * ppr_node$r[1],
      odl2 = odl2 * ppr_node$r[2],
      odl3 = odl3 * ppr_node$r[3]
    ) %>% # wazone r
    mutate(
      r = ifelse(
        is.na(r),
        (odl1 + odl2 + odl3) / suma_odl, # srednie wazone r
        r
      ),
      r = round(r, 3)
    ) %>%
    st_as_sf()
}
