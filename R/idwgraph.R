#' Network IDW Interpolation
#'
#' @param roads An sf object with LINESTRING geometry (the road network).
#' @param aprs An sf object with POINT geometry and a 'value' column (the observations).
#' @param power The IDW power parameter (default = 2).
#' @return An sf object of the network nodes with interpolated values.
#' @export
idwgraph <- function(roads, aprs, power = 2) {
  # 1. Ensure we have the necessary libraries
  # (These should be in your DESCRIPTION file)
  requireNamespace("sf", quietly = TRUE)
  requireNamespace("tidygraph", quietly = TRUE)
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("sfnetworks", quietly = TRUE) # Highly recommended for this

  # 2. Build the network graph
  # We use sfnetworks as it's the most robust way to snap points to lines
  net <- sfnetworks::as_sfnetwork(roads, directed = FALSE) %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(weight = as.numeric(sf::st_length(geometry)))

  # 3. Snap your points (aprs) to the nearest network nodes
  # We find the nearest node index for each point in 'aprs'
  nearest_node_indices <- sf::st_nearest_feature(
    aprs,
    sfnetworks::st_as_sf(net, "nodes")
  )

  # 4. Prepare the node data
  nodes_df <- net %>%
    tidygraph::activate(nodes) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(nodeID = dplyr::row_number(), r = NA_real_)

  # Assign the values from 'aprs' to the snapped nodes
  # If multiple points snap to the same node, we take the mean
  aprs_values <- aprs %>%
    sf::st_drop_geometry() %>%
    dplyr::mutate(node_idx = nearest_node_indices) %>%
    dplyr::group_by(node_idx) %>%
    dplyr::summarise(val = mean(value, na.rm = TRUE))

  nodes_df$r[aprs_values$node_idx] <- aprs_values$val

  # 5. Prepare the edge data
  edges_df <- net %>%
    tidygraph::activate(edges) %>%
    tibble::as_tibble()

  # 6. Call the Rust function
  interpolated_values <- ridwgraph(
    node_ids = as.integer(nodes_df$nodeID),
    initial_values = as.numeric(nodes_df$r),
    edges_from = as.integer(edges_df$from),
    edges_to = as.integer(edges_df$to),
    edges_weight = as.numeric(edges_df$weight),
    power = as.numeric(power)
  )

  # 7. Return the nodes as an sf object with the results
  wynik <- nodes_df %>%
    dplyr::mutate(r = interpolated_values) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = sf::st_crs(roads))

  return(wynik)
}
