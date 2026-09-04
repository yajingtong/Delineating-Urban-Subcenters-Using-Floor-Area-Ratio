#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(sf)
  library(stringr)
})

# Parallel version of 01_load_crop_building.R using:
#   1. city centers and city names from the Zipf boundary file,
#   2. Zipf-scaled center-to-periphery buffers as city boundaries,
#   3. EU building tiles in data/building/eu,
#   4. per-city parallel processing.
#
# The script keeps buildings whose polygons intersect the Zipf boundary. It does
# not cut building polygons at the boundary, so footprint/height/floor-area
# attributes remain attached to the original building polygon.

args_file <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args_file[grepl(file_arg, args_file)])
script_dir <- if (length(script_path) > 0) dirname(normalizePath(script_path[1])) else getwd()
root_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
if (!dir.exists(file.path(root_dir, "data"))) root_dir <- getwd()

building_root <- file.path(root_dir, "data", "building", "eu")
boundary_gpkg <- file.path(
  root_dir,
  "outputs",
  "zipf_scaled_city_periphery_distance",
  "zipf_scaled_city_boundary_buffers_from_paris_99p56km.gpkg"
)
boundary_layer <- "zipf_scaled_city_boundary_buffers"

output_dir <- file.path(root_dir, "outputs", "eu_buildings_zipf_boundaries")
output_city_dir <- file.path(output_dir, "gpkg_by_city")
output_summary_csv <- file.path(output_dir, "eu_buildings_zipf_boundaries_summary.csv")
output_tile_index_csv <- file.path(output_dir, "eu_building_tile_index_used.csv")
output_combined_gpkg <- file.path(output_dir, "eu_buildings_zipf_boundaries_all_cities.gpkg")
output_combined_rds <- file.path(output_dir, "eu_buildings_zipf_boundaries_all_cities.rds")

target_crs <- 3035
write_combined_outputs <- FALSE
detected_cores <- parallel::detectCores(logical = FALSE)
if (is.na(detected_cores) || detected_cores < 2) {
  detected_cores <- parallel::detectCores(logical = TRUE)
}
if (is.na(detected_cores) || detected_cores < 2) {
  detected_cores <- 2
}
workers <- max(1, min(detected_cores - 1, 6))

safe_file_stem <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^[:alnum:]]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  gsub("_+", "_", x)
}

clean_sf_field_names <- function(x) {
  geom_col <- attr(x, "sf_column")
  geom <- st_geometry(x)
  attrs <- st_drop_geometry(x)
  clean_names <- names(attrs) |>
    iconv(to = "ASCII//TRANSLIT") |>
    tolower()
  clean_names <- gsub("[^[:alnum:]_]+", "_", clean_names)
  clean_names <- gsub("^_+|_+$", "", clean_names)
  clean_names[clean_names == ""] <- "field"
  clean_names[clean_names %in% c("fid", "ogc_fid", "gpkg_fid")] <- paste0("source_", clean_names[clean_names %in% c("fid", "ogc_fid", "gpkg_fid")])
  names(attrs) <- make.unique(clean_names, sep = "_")
  st_sf(attrs, geometry = geom, crs = st_crs(x)) |>
    st_set_geometry("geometry")
}

parse_tile_bbox <- function(path) {
  stem <- tools::file_path_sans_ext(basename(path))
  parts <- str_split(stem, "_", simplify = TRUE)
  if (ncol(parts) < 5) {
    return(NULL)
  }

  nums <- suppressWarnings(as.numeric(parts[1, 2:5]))
  if (any(is.na(nums))) {
    return(NULL)
  }

  data.frame(
    file = path,
    tile_name = stem,
    xmin = nums[1],
    ymin = nums[2],
    xmax = nums[3],
    ymax = nums[4],
    stringsAsFactors = FALSE
  )
}

make_tile_index <- function(building_root) {
  building_files <- list.files(
    building_root,
    pattern = "[.]shp$",
    recursive = TRUE,
    full.names = TRUE
  )

  tile_df <- map(building_files, parse_tile_bbox) |>
    compact() |>
    bind_rows() |>
    distinct(tile_name, .keep_all = TRUE)

  if (nrow(tile_df) == 0) {
    stop("No parseable EU building shapefiles found under: ", building_root)
  }

  tile_geom <- pmap(
    list(tile_df$xmin, tile_df$ymin, tile_df$xmax, tile_df$ymax),
    \(xmin, ymin, xmax, ymax) {
      st_polygon(list(matrix(
        c(
          xmin, ymin,
          xmax, ymin,
          xmax, ymax,
          xmin, ymax,
          xmin, ymin
        ),
        ncol = 2,
        byrow = TRUE
      )))
    }
  )

  st_as_sf(tile_df, geometry = st_sfc(tile_geom, crs = 4326))
}

read_crop_one_tile <- function(tile_file, boundary_target, target_crs) {
  buildings <- st_read(tile_file, quiet = TRUE)
  if (is.na(st_crs(buildings))) {
    st_crs(buildings) <- 4326
  }
  buildings <- st_transform(buildings, target_crs)
  idx <- lengths(st_intersects(buildings, boundary_target)) > 0
  if (!any(idx)) {
    return(NULL)
  }
  buildings <- buildings[idx, ]
  clean_sf_field_names(buildings)
}

process_city <- function(i, boundaries, tile_index) {
  city_row <- boundaries[i, ]
  city_name <- as.character(city_row$city[1])
  city_stem <- safe_file_stem(city_name)

  message("Processing ", city_name)

  boundary_wgs84 <- st_transform(city_row, 4326)
  matching_tiles <- tile_index[lengths(st_intersects(tile_index, boundary_wgs84)) > 0, ]

  if (nrow(matching_tiles) == 0) {
    warning("No  building tiles overlap Zipf boundary for ", city_name)
    return(list(
      buildings = NULL,
      summary = data.frame(
        city = city_name,
        city_id = as.character(city_row$city_id[1]),
        boundary_radius_km = as.numeric(city_row$boundary_radius_km[1]),
        n_tiles = 0,
        n_buildings = 0,
        footprint_area_m2 = 0,
        floor_area_m2 = 0,
        volume_m3 = 0,
        gpkg_file = NA_character_
      )
    ))
  }

  boundary_target <- st_transform(city_row, target_crs)
  building_list <- map(matching_tiles$file, read_crop_one_tile, boundary_target, target_crs)
  building_list <- compact(building_list)

  if (length(building_list) == 0) {
    buildings_city <- st_sf(
      city = character(),
      geometry = st_sfc(crs = target_crs)
    )
  } else {
    buildings_city <- bind_rows(building_list)
  }

  if (nrow(buildings_city) > 0) {
    buildings_city <- buildings_city |>
      mutate(
        city = city_name,
        city_id = as.character(city_row$city_id[1]),
        boundary_radius_km = as.numeric(city_row$boundary_radius_km[1]),
        population_2025_center = as.numeric(city_row$population_2025[1]),
        bi = as.numeric(st_area(geometry)),
        hi = as.numeric(height),
        si = pmax(1, round(hi / 3)),
        gi = bi * si,
        vi = bi * hi
      ) |>
      filter(!st_is_empty(geometry))
  }

  city_gpkg <- file.path(output_city_dir, paste0(city_stem, "_eu_buildings_zipf_boundary.gpkg"))
  if (file.exists(city_gpkg)) unlink(city_gpkg)
  st_write(
    buildings_city,
    city_gpkg,
    layer = "eu_buildings_zipf_boundary",
    layer_options = "FID=gpkg_fid",
    quiet = TRUE
  )

  list(
    buildings = buildings_city,
    summary = data.frame(
      city = city_name,
      city_id = as.character(city_row$city_id[1]),
      boundary_radius_km = as.numeric(city_row$boundary_radius_km[1]),
      n_tiles = nrow(matching_tiles),
      n_buildings = nrow(buildings_city),
      footprint_area_m2 = if (nrow(buildings_city) > 0) sum(buildings_city$bi, na.rm = TRUE) else 0,
      floor_area_m2 = if (nrow(buildings_city) > 0) sum(buildings_city$gi, na.rm = TRUE) else 0,
      volume_m3 = if (nrow(buildings_city) > 0) sum(buildings_city$vi, na.rm = TRUE) else 0,
      gpkg_file = basename(city_gpkg)
    )
  )
}

if (!dir.exists(building_root)) {
  stop("Cannot find EU building folder: ", building_root)
}
if (!file.exists(boundary_gpkg)) {
  stop("Cannot find Zipf boundary GeoPackage: ", boundary_gpkg)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_city_dir, recursive = TRUE, showWarnings = FALSE)

message("Indexing raw EU building tiles...")
tile_index <- make_tile_index(building_root)
write_csv(st_drop_geometry(tile_index), output_tile_index_csv)
message("Indexed ", nrow(tile_index), " raw EU building tiles.")

boundaries <- st_read(boundary_gpkg, layer = boundary_layer, quiet = TRUE) |>
  st_make_valid() |>
  mutate(
    city = coalesce(city_label, center_name, Name),
    city = if_else(city == "Marseille - Aix-en-Provence", "Marseille", city),
    population_2025 = as.numeric(population_2025),
    boundary_radius_km = as.numeric(boundary_radius_km)
  ) |>
  arrange(desc(population_2025), city)

message("Processing ", nrow(boundaries), " cities with ", workers, " parallel workers.")

if (.Platform$OS.type != "windows" && workers > 1) {
  results <- parallel::mclapply(
    seq_len(nrow(boundaries)),
    process_city,
    boundaries = boundaries,
    tile_index = tile_index,
    mc.cores = workers
  )
} else {
  results <- lapply(
    seq_len(nrow(boundaries)),
    process_city,
    boundaries = boundaries,
    tile_index = tile_index
  )
}

summary_df <- bind_rows(map(results, "summary")) |>
  arrange(desc(n_buildings), city)
write_csv(summary_df, output_summary_csv)

if (write_combined_outputs) {
  combined_buildings <- map(results, "buildings") |>
    compact() |>
    bind_rows()

  if (file.exists(output_combined_gpkg)) unlink(output_combined_gpkg)
  st_write(
    combined_buildings,
    output_combined_gpkg,
    layer = "eu_buildings_zipf_boundaries_all_cities",
    layer_options = "FID=gpkg_fid",
    quiet = TRUE
  )
  saveRDS(combined_buildings, output_combined_rds)
}

message("Saved per-city building GeoPackages to: ", output_city_dir)
message("Saved summary to: ", output_summary_csv)
message("Saved tile index to: ", output_tile_index_csv)
if (!write_combined_outputs) {
  message("Skipped combined output. Set write_combined_outputs <- TRUE to enable it.")
}
print(summary_df, row.names = FALSE)
