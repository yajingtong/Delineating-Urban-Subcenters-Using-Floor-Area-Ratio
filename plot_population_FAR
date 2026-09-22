#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(sf)
})

# Plot GHSL population per 500 m cell against observed FAR for all 41 cities.

args_file <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args_file[grepl("--file=", args_file)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()
root_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root_dir, "outputs"))) root_dir <- getwd()

input_gpkg <- file.path(
  root_dir, "outputs", "far_residual_lisa_queen_41_footprint_cut",
  "far_residual_lisa_queen_41_with_ghsl_population_100m.gpkg"
)
input_layer <- "far_residual_lisa_queen_41_with_ghsl_population_100m"
output_dir <- file.path(
  root_dir, "outputs", "far_residual_lisa_queen_41_footprint_cut",
  "population_per_cell_vs_far_all_urban_areas"
)
output_png <- file.path(output_dir, "population_per_500m_cell_vs_far_all_41_urban_areas_point_plot.png")
output_pdf <- file.path(output_dir, "population_per_500m_cell_vs_far_all_41_urban_areas_point_plot.pdf")
summary_csv <- file.path(output_dir, "population_per_500m_cell_vs_far_city_statistics.csv")

if (!file.exists(input_gpkg)) stop("Missing GHSL population + FAR GeoPackage: ", input_gpkg)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cells <- st_read(input_gpkg, layer = input_layer, quiet = TRUE) |>
  st_drop_geometry() |>
  transmute(
    city = if_else(city == "Marseille - Aix-en-Provence", "Marseille", as.character(city)),
    population_2025 = as.numeric(population_2025),
    population_per_cell = as.numeric(population_500m),
    FAR = as.numeric(FAR)
  ) |>
  filter(
    !is.na(city), city != "",
    is.finite(population_per_cell), population_per_cell >= 0,
    is.finite(FAR), FAR >= 0
  )

city_order <- cells |>
  group_by(city) |>
  summarise(population_2025 = first(population_2025), .groups = "drop") |>
  arrange(desc(population_2025), city) |>
  pull(city)

if (length(city_order) != 41L) stop("Expected 41 urban areas; found ", length(city_order), ".")
cells <- cells |>
  mutate(city = factor(city, levels = city_order))

summary_df <- cells |>
  group_by(city) |>
  summarise(
    cells_plotted = n(),
    cells_with_population = sum(population_per_cell > 0),
    median_population_per_cell = median(population_per_cell),
    mean_population_per_cell = mean(population_per_cell),
    median_far = median(FAR),
    mean_far = mean(FAR),
    spearman_population_far = if_else(
      n_distinct(population_per_cell) > 1 & n_distinct(FAR) > 1,
      cor(population_per_cell, FAR, method = "spearman"),
      NA_real_
    ),
    .groups = "drop"
  ) |>
  mutate(city = as.character(city)) |>
  arrange(match(city, city_order))

write_csv(summary_df, summary_csv)

p <- ggplot(cells, aes(x = population_per_cell, y = FAR)) +
  geom_point(
    shape = 16,
    size = 0.22,
    alpha = 0.13,
    color = "#176B9A"
  ) +
  facet_wrap(vars(city), ncol = 6, drop = FALSE) +
  scale_x_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 1),
    breaks = c(0, 10, 100, 1000, 10000),
    labels = label_number(big.mark = ","),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10, sigma = 0.02),
    breaks = c(0, 0.05, 0.2, 1, 5),
    labels = label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  labs(
    x = "GHSL population per 500 m cell (pseudo-log scale)",
    y = "Observed FAR (pseudo-log scale)",
    caption = paste(
      "Cells with finite non-negative population and FAR are shown.",
      "Each point represents one 500 m cell; panels share axes."
    )
  ) +
  theme_minimal(base_family = "sans", base_size = 8.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E2E7EA", linewidth = 0.28),
    strip.background = element_rect(fill = "#EAF0F3", color = NA),
    strip.text = element_text(face = "bold", color = "#243A48", size = 8.2),
    axis.text = element_text(color = "#5A6D78", size = 6.7),
    axis.title = element_text(color = "#263B49", size = 9.5),
    legend.position = "none",
    plot.caption = element_text(color = "#667985", size = 8, hjust = 0),
    panel.spacing = grid::unit(0.65, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(output_png, p, width = 24, height = 27, units = "cm", dpi = 400, bg = "white")
ggsave(output_pdf, p, width = 24, height = 27, units = "cm", bg = "white")

message("Saved PNG: ", output_png)
message("Saved PDF: ", output_pdf)
message("Saved city statistics: ", summary_csv)
message("Urban areas plotted: ", n_distinct(cells$city))
message("Cells plotted: ", nrow(cells))
