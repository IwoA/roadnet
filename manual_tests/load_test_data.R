# Manual test script: Load geospatial test data
# This script loads geospatial data from R/test_data folder for manual testing

library(sf)
library(sfnetworks)
library(dplyr)
library(tidygraph)

# Define path to test data
test_data_path <- here::here("R", "test_data")

# Load the shapefile
test_data <- st_read(
  file.path(test_data_path, "cropt_cpt3161_cpt4071_cpt3890.shp")
)

# Display basic information about the loaded data
head(test_data)

# Load aprs

aprs <- read.csv2(
  file.path(test_data_path, "cropt_cpt3161_cpt4071_cpt3890.coef"),
  header = FALSE
)

# Source n_idw function
source(here::here("manual_tests", "n_idw.R"))


test <- n_idw(test_data, aprs)

# Source idwgraph function
source(here::here("R", "idwgraph.R"))

test_idwgraph <- idw_graph(test_data, aprs)

## APRy do mapy
aprs1 <- stringr::str_trim(aprs[1, ])
to_points <- function(x) {
  lon <- as.numeric(stringr::str_split_i(x, " ", 1))
  lat <- as.numeric(stringr::str_split_i(x, " ", 2))
  data.frame(lon = lon, lat = lat)
}

pprs <- purrr::map(aprs1, to_points) %>% purrr::list_rbind()
pprs <- st_as_sf(pprs, coords = c("lon", "lat"), crs = st_crs(test_data))
pprs$r <- as.numeric(t(aprs[2, ]))

## NaN sa bo niektore odcinki nie lacza sie z siecia drog

#Weryfikacja
library(tmap)
tmap_mode("view")
tm_shape(test_idwgraph) +
  tm_lines(col = "r", col.scale = tm_scale_continuous()) +
  tm_shape(test) +
  tm_dots(fill = "r", size = .3) + #tm_text("r") +
  tm_shape(pprs) +
  tm_dots(fill = "green", size = 1) +
  tm_text("r")
