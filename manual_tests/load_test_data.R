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

# Problem jest taki, że drogi tomtoma mają zduplikowane ID jeśli są dwukierunkowe, a idw_graph tego nie obsługuje.
# Dlatego trzeba policzyć sumy natężeń a po przeskalowaniu rozdzieliś je z powrotem wg proporcji natężeń.
test_data1 <- test_data %>%
  group_by(ID) %>%
  mutate(MEANDAY = sum(MEANDAY), WEEK = sum(WEEK)) |>
  distinct(ID, .keep_all = TRUE) |>
  ungroup()

tictoc::tic("idw_graph")
test <- idw_graph(test_data1, aprs)
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


## Weryfikacja interpolacji idw

set.seed(123)
n_repeats <- 6
sample_size <- 20

aprs_samples <- replicate(
  n_repeats,
  aprs |> slice_sample(n = sample_size),
  simplify = FALSE
)

aprs_samples <- list(
  aprs |> filter(r %in% c(17358, 12038, 15535, 18175)),
  aprs |> filter(r %in% c(12144, 34626, 28261)),
  aprs |> filter(r %in% c(14344, 9610)),
  aprs |> filter(r %in% c(14344, 10564, 16706))
)


objective_power_cv <- function(power, roads, aprs_subsets, ref) {
  losses <- vapply(
    aprs_subsets,
    function(aprs_sub) {
      pred <- idw_graph(roads, aprs_sub, power = power)
      pred_ref <- pred |>
        st_drop_geometry() |>
        select(ID, r_pred = r) |>
        left_join(ref, join_by(ID))

      max(abs(pred_ref$r - pred_ref$r_pred), na.rm = TRUE)
    },
    numeric(1)
  )

  losses <- losses[is.finite(losses)]

  if (length(losses) == 0) {
    return(Inf)
  }

  mean(losses)
}

power_grid_coarse <- seq(0.1, 2, by = 0.4)
coarse_scores <- vapply(
  power_grid_coarse,
  function(p) objective_power_cv(p, test_data1, aprs_samples, test),
  numeric(1)
)

coarse_best <- power_grid_coarse[which.min(coarse_scores)]
power_grid_fine <- seq(
  max(0.05, coarse_best - 0.3),
  min(5, coarse_best + 0.3),
  by = 0.05
)

mirai::daemons(length(power_grid_fine))

fine_scores <- purrr::map(
  power_grid_fine,
  purrr::in_parallel(
    \(p) my_func(p, test_data1, aprs_samples, test),
    my_func = objective_power_cv,
    aprs_samples = aprs_samples,
    test = test,
    test_data1 = test_data1,
    idw_graph = idw_graph
  )
)

best_power <- power_grid_fine[which.min(fine_scores)]
best_cv_loss <- min(unlist(fine_scores), na.rm = TRUE)

# Jedna próba do dalszej wizualizacji i porównania
aprs_sample <- aprs |> slice_sample(n = sample_size)
test_optim <- idw_graph(test_data, aprs, power = best_power)

test_difr <- test |> mutate(r_difr = round((r - test_optim$r) / r * 100, 3))
max_abs_difr <- max(abs(test_difr$r_difr), na.rm = TRUE)

tune_best_losses <- vapply(
  aprs_samples,
  function(aprs_sub) {
    pred <- idw_graph(test_data1, aprs_sub, power = best_power)
    pred_ref <- pred |>
      st_drop_geometry() |>
      select(ID, r_pred = r) |>
      left_join(test, join_by(ID))

    max(abs(pred_ref$r - pred_ref$r_pred), na.rm = TRUE)
  },
  numeric(1)
)

tibble(
  best_power = best_power,
  cv_mean_loss = best_cv_loss,
  holdout_loss = max_abs_difr,
  cv_median_loss = median(tune_best_losses, na.rm = TRUE),
  cv_p90_loss = unname(quantile(tune_best_losses, probs = 0.9, na.rm = TRUE))
)

tm_shape(test_difr) +
  tm_lines(
    col = "r_difr",
    col.scale = tm_scale_continuous(values = "viridis")
  ) +
  tm_shape(test_difr) +
  tm_dots(
    fill = "r_difr",
    size = .3,
    fill.scale = tm_scale_continuous(values = "viridis")
  ) + #tm_text("r") +
  tm_shape(aprs) +
  tm_dots(fill = "green", size = 1) +
  tm_text("r") +
  tm_shape(aprs_sample) +
  tm_dots(fill = "red", size = 1) +
  tm_text("r")

test_aprs <- test |> filter(r %in% aprs$r)
test_sample_aprs <- test_sample |>
  filter(ID %in% test_aprs$ID) |>
  select(ID, r) |>
  left_join(
    test_aprs |> st_drop_geometry() |> select(ID, r),
    by = "ID",
    suffix = c("_sample", "_test")
  ) |>
  mutate(difr = (r_sample - r_test) / r_test * 100)


tibble(test = test_aprs$r, test_sample = test_sample_aprs$r_sample) %>%
  filter(test != test_sample) %>%
  mutate(
    test = test - mean(test),
    test_sample = test_sample - mean(test_sample)
  ) %>%
  ggplot() +
  # geom_qq(aes(sample = test - test_sample)) +
  # geom_qq_line(aes(sample = test - test_sample)) +
  geom_point(aes(x = test, y = test_sample)) +
  geom_smooth(aes(x = test, y = test_sample), method = "lm") +
  labs(title = "Porównanie wartości r dla pełnych danych i próby") +
  theme_minimal()


tm_shape(test_sample_aprs) +
  tm_dots(
    fill = "difr",
    size = 1,
    fill.scale = tm_scale_continuous()
  )
