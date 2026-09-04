#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(sf)
  library(tibble)
})

# Compute cell-level building metrics with explicit footprint cutting:
#   1. Use the existing 500 m city cells and verified buffer_area denominator.
#   2. Cut building footprints by cell boundaries.
#   3. Allocate B, G, and V by the share of each building footprint inside each cell.
#   4. Recalculate FAR = G / buffer_area.
#
# This avoids assigning the full building to every touched cell while preserving
# the previously computed buffer_area = area(cell intersect union(buffer(buildings, 50 m))).

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_file[grepl(file_arg, args_file)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()
root_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root_dir, "outputs"))) root_dir <- getwd()

input_cells_gpkg <- file.path(
  root_dir,
  "outputs/all_city_cells_zipf_boundaries_500m_buffer",
  "all_city_cells_zipf_boundaries_500m_buffer.gpkg"
)
input_cells_layer <- "all_city_cells_zipf_boundaries_500m_buffer"

building_summary_csv <- file.path(
  root_dir,
  "outputs/eu_buildings_zipf_boundaries",
  "eu_buildings_zipf_boundaries_summary.csv"
)
building_gpkg_dir <- file.path(
  root_dir,
  "outputs/eu_buildings_zipf_boundaries",
  "gpkg_by_city"
)

output_dir <- file.path(root_dir, "outputs", "all_city_cells_zipf_boundaries_500m_buffer_footprint_cut_reused_buffer_area")
output_city_dir <- file.path(output_dir, "gpkg_by_city")
output_gpkg <- file.path(output_dir, "all_city_cells_zipf_boundaries_500m_buffer_footprint_cut_reused_buffer_area.gpkg")
output_layer <- "all_city_cells_footprint_cut_reused_buffer_area"
summary_csv <- file.path(output_dir, "all_city_cells_footprint_cut_reused_buffer_area_summary.csv")
validation_csv <- file.path(output_dir, "all_city_cells_footprint_cut_reused_buffer_area_validation.csv")

target_crs <- 3035

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_city_dir, recursive = TRUE, showWarnings = FALSE)

safe_file_stem <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  gsub("_+", "_", x)
}

allocate_building_attributes_by_cut <- function(grid_city, b_city) {
  if (nrow(grid_city) == 0 || nrow(b_city) == 0) {
    return(tibble(id = integer(), G = numeric(), V = numeric(), H = numeric(), B = numeric(), N = integer()))
  }

  b_for_allocation <- b_city |>
    mutate(b_row = row_number())
  b_for_allocation$building_area_m2 <- as.numeric(st_area(b_for_allocation))

  b_for_allocation <- b_for_allocation |>
    filter(is.finite(building_area_m2), building_area_m2 > 0) |>
    select(b_row, gi, vi, hi, bi, building_area_m2)

  parts <- tryCatch(
    suppressWarnings(st_intersection(b_for_allocation, grid_city |> select(id))),
    error = function(e) {
      warning("Footprint-cell cut failed; retrying with valid geometries: ", conditionMessage(e))
      suppressWarnings(st_intersection(
        st_make_valid(b_for_allocation),
        st_make_valid(grid_city |> select(id))
      ))
    }
  )

  if (nrow(parts) == 0) {
    return(tibble(id = integer(), G = numeric(), V = numeric(), H = numeric(), B = numeric(), N = integer()))
  }

  parts$part_area_m2 <- as.numeric(st_area(parts))

  parts |>
    mutate(
      allocation_share = if_else(
        is.finite(building_area_m2) & building_area_m2 > 0,
        pmin(part_area_m2 / building_area_m2, 1),
        0
      ),
      G_part = gi * allocation_share,
      V_part = vi * allocation_share,
      B_part = part_area_m2
    ) |>
    st_drop_geometry() |>
    group_by(id) |>
    summarise(
      G = sum(G_part, na.rm = TRUE),
      V = sum(V_part, na.rm = TRUE),
      H = if_else(
        sum(is.finite(hi) & is.finite(part_area_m2) & part_area_m2 > 0) > 0,
        weighted.mean(hi, part_area_m2, na.rm = TRUE),
        NA_real_
      ),
      B = sum(B_part, na.rm = TRUE),
      N = n_distinct(b_row),
      .groups = "drop"
    )
}

process_city <- function(city_name, cells, building_summary) {
  message("Processing city: ", city_name)

  grid_city <- cells |>
    filter(city == city_name)

  gpkg_file <- building_summary |>
    filter(city == city_name) |>
    pull(gpkg_file) |>
    first()

  if (length(gpkg_file) == 0 || is.na(gpkg_file)) {
    warning("No building GPKG listed for ", city_name)
    attr_df <- tibble(id = integer(), G = numeric(), V = numeric(), H = numeric(), B = numeric(), N = integer())
  } else {
    b_path <- file.path(building_gpkg_dir, gpkg_file)
    if (!file.exists(b_path)) {
      warning("No building GPKG found for ", city_name, ": ", b_path)
      attr_df <- tibble(id = integer(), G = numeric(), V = numeric(), H = numeric(), B = numeric(), N = integer())
    } else {
      b_city <- st_read(b_path, quiet = TRUE) |>
        st_transform(target_crs) |>
        filter(is.finite(gi), gi > 0)
      b_city <- b_city[!st_is_empty(b_city), ]
      attr_df <- allocate_building_attributes_by_cut(grid_city, b_city)
    }
  }

  out <- grid_city |>
    select(-any_of(c("G", "V", "H", "B", "N", "FAR", "far_method"))) |>
    left_join(attr_df, by = "id") |>
    mutate(
      G = coalesce(G, 0),
      V = coalesce(V, 0),
      B = coalesce(B, 0),
      N = coalesce(N, 0L),
      FAR = if_else(is.finite(buffer_area) & buffer_area > 0, G / buffer_area, NA_real_),
      far_method = "FAR = G / buffer_area; G,V,B allocated by footprint-cell intersection share; buffer_area reused from 50 m building-buffer union"
    )

  city_stem <- safe_file_stem(city_name)
  city_gpkg <- file.path(output_city_dir, paste0(city_stem, "_cells_footprint_cut_reused_buffer_area.gpkg"))
  if (file.exists(city_gpkg)) unlink(city_gpkg)
  st_write(out, city_gpkg, layer = output_layer, layer_options = "FID=gpkg_fid", quiet = TRUE)
  message("Saved ", nrow(out), " cells for ", city_name, " to ", city_gpkg)

  out
}

if (!file.exists(input_cells_gpkg)) stop("Cannot find cell GPKG: ", input_cells_gpkg)
if (!file.exists(building_summary_csv)) stop("Cannot find building summary: ", building_summary_csv)

cells <- st_read(input_cells_gpkg, layer = input_cells_layer, quiet = TRUE) |>
  st_transform(target_crs) |>
  mutate(city = as.character(city))

building_summary <- read_csv(building_summary_csv, show_col_types = FALSE) |>
  mutate(city = as.character(city)) |>
  filter(!is.na(gpkg_file), n_buildings > 0)

city_filter <- Sys.getenv("CITY_FILTER", unset = "")
cities <- sort(unique(cells$city))
if (nzchar(city_filter)) {
  selected_cities <- trimws(strsplit(city_filter, ",", fixed = TRUE)[[1]])
  cities <- intersect(cities, selected_cities)
  message("CITY_FILTER active: ", paste(selected_cities, collapse = ", "))
}

city_order <- cells |>
  st_drop_geometry() |>
  distinct(city, population_2025) |>
  filter(city %in% cities) |>
  arrange(desc(population_2025), city) |>
  pull(city)

message("Recomputing footprint-cut allocation for ", length(city_order), " cities.")

city_cells <- map(city_order, process_city, cells = cells, building_summary = building_summary)
all_city_cells <- st_as_sf(bind_rows(city_cells))

if (file.exists(output_gpkg)) unlink(output_gpkg)
st_write(all_city_cells, output_gpkg, layer = output_layer, layer_options = "FID=gpkg_fid", quiet = TRUE)

summary_tbl <- all_city_cells |>
  st_drop_geometry() |>
  group_by(city) |>
  summarise(
    n_cells = n(),
    n_cells_with_buffers = sum(is.finite(buffer_area) & buffer_area > 0, na.rm = TRUE),
    total_buffer_area = sum(buffer_area, na.rm = TRUE),
    total_G = sum(G, na.rm = TRUE),
    total_V = sum(V, na.rm = TRUE),
    total_B = sum(B, na.rm = TRUE),
    total_N = sum(N, na.rm = TRUE),
    mean_FAR = mean(FAR, na.rm = TRUE),
    median_FAR = median(FAR, na.rm = TRUE),
    max_FAR = max(FAR, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(desc(n_cells), city)

validation_tbl <- all_city_cells |>
  st_drop_geometry() |>
  filter(is.finite(buffer_area), buffer_area > 0) |>
  group_by(city) |>
  summarise(
    n_positive_buffer_cells = n(),
    max_abs_far_error = max(abs(FAR - G / buffer_area), na.rm = TRUE),
    max_buffer_area = max(buffer_area, na.rm = TRUE),
    max_cell_area = max(cell_area, na.rm = TRUE),
    n_buffer_area_gt_cell_area = sum(buffer_area > cell_area + 1e-6, na.rm = TRUE),
    n_footprint_gt_buffer_area = sum(is.finite(B) & B > buffer_area + 1e-6, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(city)

write_csv(summary_tbl, summary_csv)
write_csv(validation_tbl, validation_csv)

message("Saved combined GPKG to: ", output_gpkg)
message("Saved per-city GPKGs to: ", output_city_dir)
message("Saved summary to: ", summary_csv)
message("Saved validation to: ", validation_csv)
message("Cities processed: ", n_distinct(all_city_cells$city))
print(summary_tbl, n = nrow(summary_tbl))
print(validation_tbl, n = nrow(validation_tbl))
