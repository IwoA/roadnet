# Manual test script: Load geospatial test data
# This script loads geospatial data from R/test_data folder for manual testing

# Define path to test data
test_data_path <- here::here("R", "test_data")

# Load the shapefile
test_data <- sf::st_read(
  #file.path(test_data_path, "cropt_cpt3161_cpt4071_cpt3890.shp")
  file.path(test_data_path, "legnica_test.gpkg")
)

# Display basic information about the loaded data
head(test_data)

# Load aprs

aprs <- sf::st_read(
  file.path(test_data_path, "legnica_pprs.gpkg")
)

# Source idwgraph function
source(here::here("R", "idwgraph.R"))

tictoc::tic("idw_graph")
test <- idw_graph(test_data, aprs)
tictoc::toc()

# Weryfikacja
library(tmap)
tmap_mode("view")

tm_shape(test) +
  tm_lines(col = "r", col.scale = tm_scale_continuous(values = "viridis")) +
  tm_shape(test) +
  tm_dots(
    fill = "r",
    size = .3,
    fill.scale = tm_scale_continuous(values = "viridis")
  ) + #tm_text("r") +
  tm_shape(aprs) +
  tm_dots(fill = "green", size = 1) +
  tm_text("r")

# Łączenie z oryginalnymi danymi
test_joined <- test_data %>%
  left_join(test %>% st_drop_geometry() %>% select(ID, r), by = "ID") |>
  mutate()
