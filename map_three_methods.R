#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(tidyr)
})

sf_use_s2(FALSE)

# Compare absolute FAR thresholds, observed-FAR LISA, and signed power-law
# residual LISA for Paris in one common CBD-centred extent.

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_file[grepl(file_arg, args_file)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()
root_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root_dir, "outputs"))) root_dir <- getwd()

observed_gpkg <- file.path(
  root_dir,
  "outputs/observed_far_lisa_queen_41_footprint_cut",
  "observed_far_lisa_queen_41_footprint_cut.gpkg"
)
observed_layer <- "observed_far_lisa_queen_41"
residual_gpkg <- file.path(
  root_dir,
  "outputs/observed_vs_power_law_residual_lisa_41_cities",
  "power_law_residual_lisa_queen_41_cities.gpkg"
)
residual_layer <- "power_law_residual_lisa_queen_41_cities"
centers_gpkg <- file.path(root_dir, "outputs", "centers_20_plus_ghs_ucdb_29_merged.gpkg")
centers_layer <- "centers_20_plus_ghs_ucdb_29_merged"

output_dir <- file.path(root_dir, "outputs", "paris_far_threshold_observed_lisa_residual_lisa")
output_png <- file.path(output_dir, "paris_far_threshold_observed_lisa_residual_lisa_zoomed_full_extent.png")
output_pdf <- file.path(output_dir, "paris_far_threshold_observed_lisa_residual_lisa_zoomed_full_extent.pdf")
output_gpkg <- file.path(output_dir, "paris_far_threshold_observed_lisa_residual_lisa_zoomed_full_extent.gpkg")
output_csv <- file.path(output_dir, "paris_far_threshold_observed_lisa_residual_lisa_zoomed_full_extent_summary.csv")

for (path in c(observed_gpkg, residual_gpkg, centers_gpkg)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

target_crs <- st_crs(3035)
zoom_radius_m <- 10000

far_levels <- c("Low (0.4-<0.8)", "Medium (0.8-<2)", "High (>=2)")
far_colors <- c(
  "Low (0.4-<0.8)" = "#F3D37A",
  "Medium (0.8-<2)" = "#E88945",
  "High (>=2)" = "#B3202A"
)
lisa_levels <- c("High-High", "High-Low", "Low-High", "Low-Low", "Not Significant")
lisa_colors <- c(
  "High-High" = "#D73027",
  "High-Low" = "#E78AC3",
  "Low-High" = "#2C7FB8",
  "Low-Low" = "#084081",
  "Not Significant" = "#E1E1DE"
)

classify_far <- function(x) {
  case_when(
    !is.finite(x) | x < 0.4 ~ NA_character_,
    x < 0.8 ~ far_levels[[1]],
    x < 2 ~ far_levels[[2]],
    TRUE ~ far_levels[[3]]
  )
}

observed_query <- paste0("SELECT * FROM ", observed_layer, " WHERE city = 'Paris'")
observed <- st_read(observed_gpkg, query = observed_query, quiet = TRUE) |>
  st_transform(target_crs) |>
  st_make_valid() |>
  mutate(
    FAR_class = factor(classify_far(as.numeric(FAR)), levels = far_levels),
    lisa_cluster = recode(
      as.character(observed_far_lisa_cluster),
      "No valid FAR" = "No valid value"
    ),
    lisa_cluster = factor(lisa_cluster, levels = c(lisa_levels, "No valid value"))
  )

message("Reading Paris power-law residual LISA cells")
residual_query <- paste0("SELECT * FROM ", residual_layer, " WHERE city = 'Paris'")
residual <- st_read(residual_gpkg, query = residual_query, quiet = TRUE) |>
  st_transform(target_crs) |>
  st_make_valid() |>
  mutate(
    lisa_cluster = factor(
      as.character(power_law_residual_lisa_cluster),
      levels = c(lisa_levels, "No valid value")
    )
  )

message("Reading the Paris CBD")
center <- st_read(centers_gpkg, layer = centers_layer, quiet = TRUE) |>
  st_drop_geometry() |>
  mutate(
    city = as.character(coalesce(na_if(as.character(city), ""), requested_city)),
    center_lon = as.numeric(center_lon),
    center_lat = as.numeric(center_lat)
  ) |>
  filter(city == "Paris", is.finite(center_lon), is.finite(center_lat)) |>
  distinct(city, .keep_all = TRUE) |>
  st_as_sf(coords = c("center_lon", "center_lat"), crs = 4326, remove = FALSE) |>
  st_transform(target_crs)
if (nrow(center) != 1L) stop("Expected exactly one Paris center.")

center_xy <- st_coordinates(center)[1, ]
zoom_xlim <- center_xy[[1]] + c(-zoom_radius_m, zoom_radius_m)
zoom_ylim <- center_xy[[2]] + c(-zoom_radius_m, zoom_radius_m)
zoom_extent <- st_as_sfc(
  st_bbox(
    c(xmin = zoom_xlim[[1]], ymin = zoom_ylim[[1]], xmax = zoom_xlim[[2]], ymax = zoom_ylim[[2]]),
    crs = target_crs
  )
) |>
  st_sf(extent = "Paris CBD 20 km x 20 km", geometry = _)

observed_bbox <- st_bbox(observed)
full_radius_m <- max(
  abs(observed_bbox[["xmin"]] - center_xy[[1]]),
  abs(observed_bbox[["xmax"]] - center_xy[[1]]),
  abs(observed_bbox[["ymin"]] - center_xy[[2]]),
  abs(observed_bbox[["ymax"]] - center_xy[[2]])
)
full_radius_m <- ceiling((full_radius_m + 1000) / 5000) * 5000
full_xlim <- center_xy[[1]] + c(-full_radius_m, full_radius_m)
full_ylim <- center_xy[[2]] + c(-full_radius_m, full_radius_m)
full_extent <- st_as_sfc(
  st_bbox(
    c(xmin = full_xlim[[1]], ymin = full_ylim[[1]], xmax = full_xlim[[2]], ymax = full_ylim[[2]]),
    crs = target_crs
  )
) |>
  st_sf(extent = "Full Paris study extent", geometry = _)

observed_zoom <- observed[lengths(st_intersects(observed, zoom_extent)) > 0, ]
residual_zoom <- residual[lengths(st_intersects(residual, zoom_extent)) > 0, ]
if (nrow(observed_zoom) == 0L || nrow(residual_zoom) == 0L) stop("No cells intersect the Paris zoom extent.")

summarise_extent <- function(observed_data, residual_data, extent_name) {
  bind_rows(
    observed_data |>
      st_drop_geometry() |>
      mutate(category = as.character(FAR_class)) |>
      filter(!is.na(category)) |>
      count(category, name = "cells") |>
      mutate(panel = "FAR threshold"),
    observed_data |>
      st_drop_geometry() |>
      mutate(category = as.character(lisa_cluster)) |>
      filter(category %in% lisa_levels) |>
      count(category, name = "cells") |>
      mutate(panel = "Observed FAR LISA"),
    residual_data |>
      st_drop_geometry() |>
      mutate(category = as.character(lisa_cluster)) |>
      filter(category %in% lisa_levels) |>
      count(category, name = "cells") |>
      mutate(panel = "Power-law residual LISA")
  ) |>
    mutate(extent = extent_name) |>
    group_by(extent, panel) |>
    mutate(percent = 100 * cells / sum(cells)) |>
    ungroup() |>
    select(extent, panel, category, cells, percent)
}

summary_df <- bind_rows(
  summarise_extent(observed_zoom, residual_zoom, "Zoomed 20 km x 20 km"),
  summarise_extent(observed, residual, "Full extent")
)
write.csv(summary_df, output_csv, row.names = FALSE, fileEncoding = "UTF-8")

base_map <- function(title, background_data, x_limits, y_limits) {
  ggplot() +
    geom_sf(data = background_data, fill = "#F1F0EC", color = NA) +
    coord_sf(xlim = x_limits, ylim = y_limits, expand = FALSE, datum = NA) +
    scale_shape_manual(values = c("CBD" = 21), name = NULL) +
    labs(title = title) +
    theme_void(base_size = 10, base_family = "serif") +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.background = element_rect(fill = "white", color = "#AFAFAF", linewidth = 0.35),
      panel.background = element_rect(fill = "#FAFAF7", color = NA),
      legend.position = "bottom",
      plot.margin = margin(3, 3, 3, 3)
    )
}

map_features <- function(x_limits, y_limits, scale_bar_m) {
  scale_x <- x_limits[[1]] + 0.065 * diff(x_limits)
  scale_y <- y_limits[[1]] + 0.065 * diff(y_limits)
  scale_tick <- 0.008 * diff(y_limits)
  center_size <- if (scale_bar_m >= 20000) 1.7 else 2.3
  text_size <- if (scale_bar_m >= 20000) 2.3 else 2.7

  list(
    geom_sf(
      data = center,
      aes(shape = "CBD"),
      fill = "white", color = "#111111", size = center_size, stroke = 0.55
    ),
    annotate(
      "segment", x = scale_x, xend = scale_x + scale_bar_m,
      y = scale_y, yend = scale_y,
      color = "#171717", linewidth = 0.65, lineend = "butt"
    ),
    annotate(
      "segment",
      x = c(scale_x, scale_x + scale_bar_m), xend = c(scale_x, scale_x + scale_bar_m),
      y = scale_y - scale_tick, yend = scale_y + scale_tick,
      color = "#171717", linewidth = 0.3
    ),
    annotate(
      "text", x = scale_x + scale_bar_m / 2, y = scale_y + 2.3 * scale_tick,
      label = paste0(scale_bar_m / 1000, " km"),
      family = "serif", fontface = "bold", size = text_size
    )
  )
}

make_threshold_panel <- function(data, title, x_limits, y_limits, scale_bar_m) {
  base_map(title, data, x_limits, y_limits) +
    geom_sf(
      data = data |> filter(!is.na(FAR_class)),
      aes(fill = FAR_class), color = NA
    ) +
    map_features(x_limits, y_limits, scale_bar_m) +
    scale_fill_manual(
      values = far_colors,
      breaks = far_levels,
      drop = FALSE,
      name = "FAR category"
    ) +
    guides(
      fill = guide_legend(order = 1, nrow = 1),
      shape = guide_legend(order = 2, override.aes = list(fill = "white", size = 3))
    )
}

make_lisa_panel <- function(background_data, data, title, x_limits, y_limits, scale_bar_m) {
  base_map(title, background_data, x_limits, y_limits) +
    geom_sf(
      data = data |> filter(as.character(lisa_cluster) %in% lisa_levels),
      aes(fill = lisa_cluster), color = NA
    ) +
    map_features(x_limits, y_limits, scale_bar_m) +
    scale_fill_manual(
      values = lisa_colors,
      breaks = lisa_levels,
      limits = lisa_levels,
      drop = FALSE,
      name = "LISA cluster"
    ) +
    guides(
      fill = guide_legend(order = 1, nrow = 1),
      shape = guide_legend(order = 2, override.aes = list(fill = "white", size = 3))
    )
}

zoom_threshold_panel <- make_threshold_panel(
  observed_zoom, "FAR threshold", zoom_xlim, zoom_ylim, 5000
)
zoom_observed_panel <- make_lisa_panel(
  observed_zoom, observed_zoom, "Observed FAR LISA", zoom_xlim, zoom_ylim, 5000
)
zoom_residual_panel <- make_lisa_panel(
  observed_zoom, residual_zoom, "Power-law residual LISA", zoom_xlim, zoom_ylim, 5000
)

full_threshold_panel <- make_threshold_panel(
  observed, "FAR threshold", full_xlim, full_ylim, 20000
)
full_observed_panel <- make_lisa_panel(
  observed, observed, "Observed FAR LISA", full_xlim, full_ylim, 20000
)
full_residual_panel <- make_lisa_panel(
  observed, residual, "Power-law residual LISA", full_xlim, full_ylim, 20000
)

threshold_legend <- cowplot::get_legend(
  zoom_threshold_panel +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", family = "serif", size = 9.5),
      legend.text = element_text(family = "serif", size = 8.8)
    )
)
lisa_legend <- cowplot::get_legend(
  zoom_residual_panel +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold", family = "serif", size = 9.5),
      legend.text = element_text(family = "serif", size = 8.8)
    )
)
legend_row <- cowplot::plot_grid(threshold_legend, lisa_legend, nrow = 1, rel_widths = c(0.9, 1.5))

zoom_row <- wrap_plots(
  lapply(
    list(zoom_threshold_panel, zoom_observed_panel, zoom_residual_panel),
    function(plot) plot + theme(legend.position = "none")
  ),
  nrow = 1
) +
  plot_annotation(
    title = "Zoomed CBD extent (20 km x 20 km)",
    theme = theme(plot.title = element_text(family = "serif", face = "bold", size = 12, hjust = 0))
  )

full_row <- wrap_plots(
  lapply(
    list(full_threshold_panel, full_observed_panel, full_residual_panel),
    function(plot) plot + theme(legend.position = "none")
  ),
  nrow = 1
) +
  plot_annotation(
    title = "Full Paris study extent",
    theme = theme(plot.title = element_text(family = "serif", face = "bold", size = 12, hjust = 0))
  )

final_plot <- wrap_plots(
  list(
    wrap_elements(full = zoom_row),
    wrap_elements(full = full_row),
    wrap_elements(full = legend_row)
  ),
  ncol = 1,
  heights = c(7, 7, 0.8)
) +
  plot_annotation(
    title = "Paris: FAR thresholds and spatial clusters",
    theme = theme(
      plot.title = element_text(family = "serif", face = "bold", size = 17, hjust = 0.5)
    )
  )

unlink(c(output_gpkg, paste0(output_gpkg, "-wal"), paste0(output_gpkg, "-shm")), force = TRUE)
st_write(
  observed |> filter(!is.na(FAR_class)),
  output_gpkg, layer = "far_threshold_cells_full", quiet = TRUE
)
st_write(
  observed |> filter(as.character(lisa_cluster) != "No valid value"),
  output_gpkg, layer = "observed_far_lisa_full", append = TRUE, quiet = TRUE
)
st_write(
  residual |> filter(as.character(lisa_cluster) != "No valid value"),
  output_gpkg, layer = "power_law_residual_lisa_full", append = TRUE, quiet = TRUE
)
st_write(center, output_gpkg, layer = "paris_cbd", append = TRUE, quiet = TRUE)
st_write(zoom_extent, output_gpkg, layer = "map_extent_zoomed", append = TRUE, quiet = TRUE)
st_write(full_extent, output_gpkg, layer = "map_extent_full", append = TRUE, quiet = TRUE)

ggsave(output_png, final_plot, width = 16, height = 12.5, dpi = 400, bg = "white")
ggsave(
  output_pdf, final_plot,
  width = 16, height = 12.5,
  device = grDevices::pdf, useDingbats = FALSE, bg = "white"
)
