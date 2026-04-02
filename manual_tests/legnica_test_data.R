library(tidyverse)
library(sf)

pprs <- readr::read_csv2(
  paste0(
    "C:\\Users\\iwo.augustynski\\Documents\\CUPT\\Data\\ITS\\Legnica\\",
    "Sumaryczne natężenie ruchu na monitorowanych skrzyżowaniach.csv"
  ),
  locale = readr::locale(encoding = "ISO-8859-2")
)

pprs <- separate(pprs, col = "Lokalizacja", into = c("y", "x"), sep = ", ") |>
  st_as_sf(coords = c("x", "y"), crs = 4326)

colnames(pprs) <- c("Lokalizacja", "r", "l_poj_25", "geometry")

pprs <- pprs |> mutate(r = round(r / 365))

doln <- read_sf(
  "C:\\Users\\iwo.augustynski\\Documents\\CUPT\\Data\\Wyniki\\tomtom_2024_dolnoslaskie.gpkg"
)
legnica <- read_sf(
  "C:\\\\Users\\iwo.augustynski\\Documents\\CUPT\\Data\\woj_shp\\Granice\\A03_Granice_gmin.shp"
) |>
  filter(JPT_NAZWA_ == "Legnica") |>
  select(JPT_NAZWA_) |>
  st_transform(crs = st_crs(doln))

test_legnica <- st_filter(doln, legnica)
test_legnica <- test_legnica |>
  rename(geometry = geom) |>
  mutate(dl_segm = st_length(geometry))

write_sf(test_legnica, "./R/test_data/legnica_test.gpkg")
write_sf(pprs, "./R/test_data/legnica_pprs.gpkg")
