#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cowplot)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(stringr)
  library(terra)
  library(tidyr)
})

sf_use_s2(FALSE)

# Map every merged HH subcenter patch in all 41 cities using four mutually
# exclusive functional categories. Patches outside France are removed and
# border-crossing patches are clipped to the national boundary. Buildings are
# rasterized to 200 m and cached because individual city GeoPackages can
# contain hundreds of thousands of polygons.

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_file[grepl(file_arg, args_file)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()
root_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root_dir, "outputs"))) root_dir <- getwd()

target_crs <- st_crs(3035)
building_resolution_m <- 200

patch_gpkg <- file.path(
  root_dir,
  "outputs/far_power_law_positive_residual_hh_patches_41",
  "subcenter_3km_office_occupancy_table",
  "retained_subcenters_office_occupancy_3km_41_cities.gpkg"
)
patch_layer <- "all_patches_3km_classified"
insee_gpkg <- file.path(root_dir, "outputs", "insee_41_urban_area_boundaries_aav2020_2026.gpkg")
insee_layer <- "insee_41_urban_area_boundaries"
centers_gpkg <- file.path(root_dir, "outputs", "centers_20_plus_ghs_ucdb_29_merged.gpkg")
centers_layer <- "centers_20_plus_ghs_ucdb_29_merged"
france_boundary_path <- file.path(root_dir, "data", "boundary", "FRA_adm3.shp")
building_dir <- file.path(root_dir, "outputs", "eu_buildings_zipf_boundaries", "gpkg_by_city")
building_layer <- "eu_buildings_zipf_boundary"

output_dir <- file.path(
  root_dir,
  "outputs/far_power_law_positive_residual_hh_patches_41",
  "all_41_city_subcenter_function_maps"
)
cache_dir <- file.path(output_dir, paste0("building_grid_cache_all_subcenters_", building_resolution_m, "m"))
output_png <- file.path(output_dir, "all_41_cities_all_subcenters_categories_clipped_to_france.png")
output_pdf <- file.path(output_dir, "all_41_cities_all_subcenters_categories_clipped_to_france.pdf")
output_gpkg <- file.path(output_dir, "all_41_cities_all_subcenters_categories_clipped_to_france.gpkg")
output_csv <- file.path(output_dir, "all_41_cities_all_subcenters_category_summary_clipped_to_france.csv")
output_audit_csv <- file.path(output_dir, "all_subcenters_france_clip_audit_by_city.csv")

for (path in c(patch_gpkg, insee_gpkg, centers_gpkg, france_boundary_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

normalise_city <- function(x) {
  x <- as.character(x)
  if_else(grepl("^Marseille", x), "Marseille", x)
}

safe_file_city <- function(city) {
  special <- c("Nîmes" = "N_imes", "Orléans" = "Orl_eans")
  if (city %in% names(special)) return(unname(special[[city]]))
  city |>
    stringi::stri_trans_general("Latin-ASCII") |>
    str_replace_all("[^A-Za-z0-9]+", "_") |>
    str_replace_all("^_|_$", "")
}

safe_ratio <- function(value, threshold) {
  case_when(
    is.finite(value) & is.finite(threshold) & threshold > 0 ~ value / threshold - 1,
    is.finite(value) & value > 0 & threshold == 0 ~ Inf,
    TRUE ~ -Inf
  )
}

sql_bbox_clause <- function(bb) {
  xmin <- format(unname(bb[["xmin"]]), scientific = FALSE, trim = TRUE)
  xmax <- format(unname(bb[["xmax"]]), scientific = FALSE, trim = TRUE)
  ymin <- format(unname(bb[["ymin"]]), scientific = FALSE, trim = TRUE)
  ymax <- format(unname(bb[["ymax"]]), scientific = FALSE, trim = TRUE)
  paste0(
    "r.minx <= ", xmax, " AND r.maxx >= ", xmin,
    " AND r.miny <= ", ymax, " AND r.maxy >= ", ymin
  )
}

building_file <- function(city) {
  file.path(building_dir, paste0(safe_file_city(city), "_eu_buildings_zipf_boundary.gpkg"))
}

read_building_grid <- function(city, bb) {
  cache_file <- file.path(cache_dir, paste0(safe_file_city(city), "_building_grid.rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))

  path <- building_file(city)
  if (!file.exists(path)) stop("Missing building GeoPackage for ", city, ": ", path)
  query <- paste0(
    "SELECT b.gpkg_fid, b.geom FROM ", building_layer, " AS b ",
    "JOIN rtree_", building_layer, "_geom AS r ON b.gpkg_fid = r.id ",
    "WHERE ", sql_bbox_clause(bb)
  )
  buildings <- suppressWarnings(st_read(path, query = query, quiet = TRUE))

  if (nrow(buildings) == 0) {
    result <- list(grid = data.frame(x = numeric(), y = numeric()), n_buildings = 0L)
    saveRDS(result, cache_file)
    return(result)
  }

  buildings <- st_transform(buildings, target_crs)
  template <- rast(
    xmin = bb[["xmin"]], xmax = bb[["xmax"]],
    ymin = bb[["ymin"]], ymax = bb[["ymax"]],
    resolution = building_resolution_m,
    crs = target_crs$wkt
  )
  building_raster <- rasterize(
    vect(buildings), template,
    field = 1, background = NA, touches = TRUE
  )
  result <- list(
    grid = as.data.frame(building_raster, xy = TRUE, na.rm = TRUE) |> select(x, y),
    n_buildings = nrow(buildings)
  )
  saveRDS(result, cache_file)
  result
}

message("Reading and classifying every subcenter patch")
patches_all <- st_read(patch_gpkg, layer = patch_layer, quiet = TRUE) |>
  st_transform(target_crs) |>
  st_make_valid() |>
  mutate(
    city = normalise_city(city),
    shop_exceedance = safe_ratio(patch_shops_per_cell, city_median_patch_shops_per_cell),
    office_exceedance = safe_ratio(patch_mean_offices_per_cell, city_median_patch_offices_per_cell),
    occupancy_exceedance = safe_ratio(patch_average_occupancy, city_median_patch_occupancy),
    strongest_single_function = case_when(
      office_exceedance >= shop_exceedance & office_exceedance >= occupancy_exceedance ~ "Office",
      shop_exceedance >= office_exceedance & shop_exceedance >= occupancy_exceedance ~ "Commercial",
      TRUE ~ "Residential"
    ),
    map_category = if_else(
      above_median_all_three,
      "Mixed",
      strongest_single_function
    ),
    map_category = factor(map_category, levels = c("Mixed", "Commercial", "Residential", "Office")),
    map_classification_rule = paste(
      "Mixed = above all three all-patch city medians; otherwise category",
      "with largest normalized value relative to shops, offices, occupancy median"
    )
  )

city_order <- patches_all |>
  st_drop_geometry() |>
  distinct(city, population_2025) |>
  arrange(desc(population_2025), city) |>
  pull(city)

if (length(city_order) != 41L) stop("Expected 41 cities; found ", length(city_order), ".")

message("Building dissolved France boundary and clipping subcenters")
france_boundary <- st_read(france_boundary_path, quiet = TRUE) |>
  st_transform(target_crs) |>
  st_make_valid() |>
  summarise(geometry = st_union(geometry))

intersects_france <- lengths(st_intersects(patches_all, france_boundary)) > 0
fully_inside_france <- rep(FALSE, nrow(patches_all))
fully_inside_france[intersects_france] <- lengths(
  st_within(patches_all[intersects_france, ], france_boundary)
) > 0

patches_all <- patches_all |>
  mutate(
    france_clip_status = case_when(
      !intersects_france ~ "removed_outside_france",
      fully_inside_france ~ "inside_france",
      TRUE ~ "clipped_at_france_border"
    )
  )

inside_patches <- patches_all[intersects_france & fully_inside_france, ]
edge_patches <- patches_all[intersects_france & !fully_inside_france, ]
if (nrow(edge_patches) > 0) {
  edge_patches <- suppressWarnings(st_intersection(edge_patches, france_boundary))
}
patches <- st_as_sf(bind_rows(inside_patches, edge_patches))

clip_audit <- patches_all |>
  st_drop_geometry() |>
  transmute(
    city,
    source_subcenters = 1L,
    removed_outside_france = as.integer(!intersects_france),
    clipped_at_france_border = as.integer(intersects_france & !fully_inside_france)
  ) |>
  group_by(city) |>
  summarise(
    source_subcenters = sum(source_subcenters),
    removed_outside_france = sum(removed_outside_france),
    clipped_at_france_border = sum(clipped_at_france_border),
    mapped_subcenters = source_subcenters - removed_outside_france,
    .groups = "drop"
  ) |>
  arrange(match(city, city_order))

if (length(unique(patches$city)) != 41L) stop("France clipping removed every patch from at least one city.")
message(
  "Mapping ", nrow(patches), " subcenters across 41 cities; removed ",
  sum(!intersects_france), " patches outside France and clipped ",
  sum(intersects_france & !fully_inside_france), " border patches"
)

insee_boundaries <- st_read(insee_gpkg, layer = insee_layer, quiet = TRUE) |>
  st_transform(target_crs) |>
  st_make_valid() |>
  mutate(city = normalise_city(city)) |>
  filter(city %in% city_order)

city_centers <- st_read(centers_gpkg, layer = centers_layer, quiet = TRUE) |>
  st_drop_geometry() |>
  mutate(
    city = normalise_city(coalesce(na_if(as.character(city), ""), requested_city)),
    center_lon = as.numeric(center_lon),
    center_lat = as.numeric(center_lat)
  ) |>
  filter(city %in% city_order, is.finite(center_lon), is.finite(center_lat)) |>
  distinct(city, .keep_all = TRUE) |>
  st_as_sf(coords = c("center_lon", "center_lat"), crs = 4326, remove = FALSE) |>
  st_transform(target_crs)

for (object_name in c("insee_boundaries", "city_centers")) {
  missing_cities <- setdiff(city_order, get(object_name)$city)
  if (length(missing_cities) > 0) {
    stop(object_name, " is missing: ", paste(missing_cities, collapse = ", "))
  }
}

missing_buildings <- city_order[!file.exists(vapply(city_order, building_file, character(1)))]
if (length(missing_buildings) > 0) {
  stop("Missing building files: ", paste(missing_buildings, collapse = ", "))
}

category_colors <- c(
  "Mixed" = "#C53B32",
  "Commercial" = "#E78AC3",
  "Residential" = "#F2C94C",
  "Office" = "#2C7FB8"
)

city_bbox <- function(city_name) {
  boundary_city <- insee_boundaries |> filter(city == city_name)
  patch_city <- patches |> filter(city == city_name)
  bounds <- rbind(
    st_bbox(boundary_city),
    if (nrow(patch_city) > 0) st_bbox(patch_city) else st_bbox(boundary_city)
  )
  bb <- c(
    xmin = min(bounds[, "xmin"]), ymin = min(bounds[, "ymin"]),
    xmax = max(bounds[, "xmax"]), ymax = max(bounds[, "ymax"])
  )
  padding <- max(750, 0.03 * max(bb[["xmax"]] - bb[["xmin"]], bb[["ymax"]] - bb[["ymin"]]))
  c(
    xmin = bb[["xmin"]] - padding, ymin = bb[["ymin"]] - padding,
    xmax = bb[["xmax"]] + padding, ymax = bb[["ymax"]] + padding
  )
}

make_city_plot <- function(city_name, building_grid) {
  boundary_city <- insee_boundaries |> filter(city == city_name)
  patch_city <- patches |> filter(city == city_name)
  center_city <- city_centers |> filter(city == city_name)
  bb <- city_bbox(city_name)
  width <- bb[["xmax"]] - bb[["xmin"]]
  height <- bb[["ymax"]] - bb[["ymin"]]
  scale_x <- bb[["xmin"]] + 0.05 * width
  scale_y <- bb[["ymin"]] + 0.05 * height
  scale_tick <- 0.012 * height

  legend_keys <- data.frame(
    x = rep(bb[["xmin"]], 4), y = rep(bb[["ymin"]], 4),
    map_category = factor(names(category_colors), levels = names(category_colors))
  )

  ggplot() +
    geom_point(
      data = legend_keys,
      aes(x = x, y = y, fill = map_category),
      shape = 22, size = 0, alpha = 0, show.legend = TRUE
    ) +
    geom_tile(
      data = building_grid,
      aes(x = x, y = y),
      width = building_resolution_m, height = building_resolution_m,
      fill = "#BDC4C8", alpha = 0.5
    ) +
    geom_sf(
      data = boundary_city,
      aes(color = "INSEE urban-area boundary"),
      fill = NA, linewidth = 0.32
    ) +
    geom_sf(
      data = patch_city,
      aes(fill = map_category),
      color = "#252525", linewidth = 0.04, alpha = 0.96
    ) +
    geom_sf(
      data = center_city,
      aes(shape = "CBD"),
      fill = "white", color = "#111111", size = 1.25, stroke = 0.35
    ) +
    annotate(
      "segment", x = scale_x, xend = scale_x + 20000,
      y = scale_y, yend = scale_y,
      color = "#171717", linewidth = 0.42, lineend = "butt"
    ) +
    annotate(
      "segment",
      x = c(scale_x, scale_x + 20000), xend = c(scale_x, scale_x + 20000),
      y = scale_y - scale_tick, yend = scale_y + scale_tick,
      color = "#171717", linewidth = 0.25
    ) +
    annotate(
      "text", x = scale_x + 10000, y = scale_y + 2.2 * scale_tick,
      label = "20 km", family = "serif", fontface = "bold", size = 1.75
    ) +
    coord_sf(
      xlim = c(bb[["xmin"]], bb[["xmax"]]),
      ylim = c(bb[["ymin"]], bb[["ymax"]]),
      expand = FALSE, datum = NA
    ) +
    scale_fill_manual(
      values = category_colors,
      breaks = names(category_colors),
      drop = FALSE,
      name = "Subcenter category"
    ) +
    scale_color_manual(values = c("INSEE urban-area boundary" = "#161616"), name = NULL) +
    scale_shape_manual(values = c("CBD" = 21), name = "Map feature") +
    guides(
      fill = guide_legend(
        order = 1, nrow = 1,
        override.aes = list(shape = 22, size = 3.2, color = NA, alpha = 1)
      ),
      color = guide_legend(order = 2, override.aes = list(linewidth = 0.9)),
      shape = guide_legend(order = 3, override.aes = list(fill = "white", size = 2.7))
    ) +
    labs(title = city_name) +
    theme_void(base_size = 7, base_family = "serif") +
    theme(
      plot.title = element_text(size = 8.4, face = "bold", hjust = 0.5, color = "#202020", margin = margin(b = 1)),
      plot.background = element_rect(fill = "white", color = "#BDBDBD", linewidth = 0.2),
      panel.background = element_rect(fill = "#FAFAF8", color = NA),
      legend.position = "bottom",
      legend.title = element_text(size = 8.5, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(2, 2, 2, 2)
    )
}

message("Loading cached/building backgrounds and creating ", length(city_order), " panels")
plots <- vector("list", length(city_order))
building_counts <- setNames(integer(length(city_order)), city_order)

for (i in seq_along(city_order)) {
  city_name <- city_order[[i]]
  message("  [", i, "/", length(city_order), "] ", city_name)
  bb <- city_bbox(city_name)
  building_result <- read_building_grid(city_name, bb)
  building_counts[[city_name]] <- building_result$n_buildings
  plots[[i]] <- make_city_plot(city_name, building_result$grid)
  rm(building_result)
  invisible(gc(FALSE))
}

shared_legend <- cowplot::get_legend(
  plots[[1]] + theme(legend.position = "bottom", legend.box = "horizontal")
)
map_panel <- wrap_plots(
  lapply(plots, function(plot) plot + theme(legend.position = "none")),
  ncol = 7
)
combined_map <- map_panel / wrap_elements(full = shared_legend) +
  plot_layout(heights = c(16, 0.75))

summary_df <- patches |>
  st_drop_geometry() |>
  count(city, map_category, name = "subcenters") |>
  complete(
    city = city_order,
    map_category = factor(names(category_colors), levels = names(category_colors)),
    fill = list(subcenters = 0L)
  ) |>
  arrange(match(city, city_order), map_category) |>
  mutate(n_buildings_mapped = building_counts[city]) |>
  left_join(clip_audit, by = "city")

unlink(c(output_gpkg, paste0(output_gpkg, "-wal"), paste0(output_gpkg, "-shm")), force = TRUE)
st_write(patches, output_gpkg, layer = "all_subcenter_categories_clipped_to_france", quiet = TRUE)
st_write(insee_boundaries, output_gpkg, layer = "insee_urban_area_boundaries", append = TRUE, quiet = TRUE)
st_write(france_boundary, output_gpkg, layer = "france_boundary", append = TRUE, quiet = TRUE)
write.csv(summary_df, output_csv, row.names = FALSE, fileEncoding = "UTF-8")
write.csv(clip_audit, output_audit_csv, row.names = FALSE, fileEncoding = "UTF-8")

ggsave(output_png, combined_map, width = 24, height = 20, dpi = 300, bg = "white")
ggsave(
  output_pdf, combined_map,
  width = 24, height = 20,
  device = grDevices::pdf, bg = "white", useDingbats = FALSE
)

message("Saved PNG: ", output_png)
message("Saved PDF: ", output_pdf)
message("Saved mapped data: ", output_gpkg)
message("Saved summary: ", output_csv)
message("Saved France clipping audit: ", output_audit_csv)
