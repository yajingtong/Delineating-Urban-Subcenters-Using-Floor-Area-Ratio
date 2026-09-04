library(sf)
library(dplyr)
library(purrr)
library(tibble)

centers_rds <- "data/centers.rds"
buildings_path <- "outputs/all_city_M.shp"
grid_rds <- "outputs/all_city_grids_100rings_full_500m.rds"
output_gpkg <- "data/cells_explicit_attribute_allocation_500m.gpkg"
output_layer <- "all_city_cells_explicit_attribute_allocation_500m"
output_csv <- "outputs/all_city_cells_explicit_attribute_allocation_500m_summary.csv"

cell_size <- 500
buffer_dist <- 50
ring_width <- 250
max_rings <- 100
radius_max <- ring_width * max_rings
target_crs <- 32632

for (f in c(centers_rds, buildings_path)) {
  if (!file.exists(f)) {
    stop("Cannot find input file: ", f)
  }
}

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

if (file.exists(output_gpkg)) {
  file.remove(output_gpkg)
}
unlink(c(output_gpkg, paste0(output_gpkg, "-wal"), paste0(output_gpkg, "-shm")), force = TRUE)

centers <- readRDS(centers_rds) %>%
  st_as_sf() %>%
  st_transform(target_crs) %>%
  mutate(Name = as.character(Name))

read_buildings <- function(path) {
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    readRDS(path) %>% st_as_sf()
  } else {
    st_read(path, quiet = TRUE)
  }
}

all_buildings <- read_buildings(buildings_path) %>%
  st_as_sf() %>%
  st_transform(target_crs) %>%
  mutate(
    city = as.character(city),
    gi = as.numeric(gi),
    vi = as.numeric(vi),
    bi = as.numeric(bi),
    hi = as.numeric(hi)
  )

create_empty_grid <- function(city_name, center_geom) {
  message("Creating grid for ", city_name)

  city_circle <- st_buffer(center_geom, radius_max)

  grid <- st_make_grid(
    city_circle,
    cellsize = c(cell_size, cell_size),
    what = "polygons"
  )

  grid_sf <- st_sf(
    id = seq_along(grid),
    geometry = grid
  ) %>%
    st_intersection(city_circle)

  cent <- st_coordinates(st_centroid(grid_sf))
  center_xy <- st_coordinates(center_geom)[1, ]

  grid_sf %>%
    mutate(
      X = cent[, 1],
      Y = cent[, 2],
      distCBD = sqrt((X - center_xy[1])^2 + (Y - center_xy[2])^2),
      ring_id = pmax(1L, ceiling(distCBD / ring_width)),
      city = city_name
    )
}

if (file.exists(grid_rds)) {
  message("Using full 100-ring grid: ", grid_rds)
  all_city_grids <- readRDS(grid_rds) %>%
    st_as_sf() %>%
    st_transform(target_crs) %>%
    mutate(city = as.character(city))
} else {
  message("Full 100-ring grid not found; creating fallback 100-ring circular grid.")
  all_city_grids <- map_dfr(
    seq_len(nrow(centers)),
    function(i) {
      create_empty_grid(
        city_name = centers$Name[i],
        center_geom = centers[i, ]
      )
    }
  )
}

allocate_building_attributes <- function(grid_city, b_city) {
  # Use building footprints for the numerator. Buffers are only used later for
  # the denominator, so a tiny buffer sliver does not pull a full building G into
  # an otherwise empty cell.
  hits <- st_intersects(b_city, grid_city)

  if (length(hits) == 0 || sum(lengths(hits)) == 0) {
    return(tibble(
      id = integer(),
      G = numeric(),
      V = numeric(),
      H = numeric(),
      B = numeric(),
      N = integer()
    ))
  }

  hit_tbl <- tibble(
    b_row = rep(seq_along(hits), lengths(hits)),
    id = unlist(hits)
  )

  building_attr <- b_city %>%
    st_drop_geometry() %>%
    mutate(b_row = row_number())

  hit_tbl %>%
    left_join(building_attr, by = "b_row") %>%
    group_by(id) %>%
    summarise(
      G = sum(gi, na.rm = TRUE),
      V = sum(vi, na.rm = TRUE),
      H = mean(hi, na.rm = TRUE),
      B = sum(bi, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    )
}

compute_buffer_area_by_cell <- function(grid_city, b_buffered) {
  inter <- st_intersection(
    grid_city %>% select(id),
    b_buffered %>% select(b_row)
  )

  if (nrow(inter) == 0) {
    return(tibble(
      id = integer(),
      buffer_area = numeric()
    ))
  }

  inter %>%
    group_by(id) %>%
    summarise(
      geometry = st_union(geometry),
      .groups = "drop"
    ) %>%
    mutate(
      buffer_area = as.numeric(st_area(geometry))
    ) %>%
    st_drop_geometry()
}

process_city_grid_explicit <- function(city_name, all_city_grids, buildings_sf, buffer_dist) {
  message("Processing city: ", city_name)

  grid_city <- all_city_grids %>%
    filter(city == city_name)

  b_city <- buildings_sf %>%
    filter(city == city_name)

  if (nrow(grid_city) == 0) {
    return(NULL)
  }

  if (nrow(b_city) == 0) {
    return(
      grid_city %>%
        mutate(
          buffer_area = NA_real_,
          G = NA_real_,
          V = NA_real_,
          H = NA_real_,
          B = NA_real_,
          N = 0L,
          FAR = NA_real_
        )
    )
  }

  b_city <- b_city %>%
    mutate(b_row = row_number())

  b_buffered <- b_city %>%
    st_buffer(dist = buffer_dist)

  attr_df <- allocate_building_attributes(grid_city, b_city)
  area_df <- compute_buffer_area_by_cell(grid_city, b_buffered)

  stats <- full_join(area_df, attr_df, by = "id") %>%
    mutate(
      buffer_area = coalesce(buffer_area, 0),
      G = coalesce(G, 0),
      V = coalesce(V, 0),
      B = coalesce(B, 0),
      N = coalesce(N, 0L),
      FAR = if_else(buffer_area > 0, G / buffer_area, NA_real_)
    )

  grid_city %>%
    left_join(stats, by = "id")
}

cities <- sort(unique(all_city_grids$city))

results_list <- map(
  cities,
  ~ process_city_grid_explicit(
    city_name = .x,
    all_city_grids = all_city_grids,
    buildings_sf = all_buildings,
    buffer_dist = buffer_dist
  )
)

all_city_cells <- bind_rows(results_list) %>%
  st_as_sf()

st_write(all_city_cells, output_gpkg, layer = output_layer, delete_dsn = TRUE, quiet = TRUE)

summary_df <- all_city_cells %>%
  st_drop_geometry() %>%
  group_by(city) %>%
  summarise(
    n_cells = n(),
    n_cells_with_buffers = sum(!is.na(buffer_area) & buffer_area > 0),
    total_buffer_area = sum(buffer_area, na.rm = TRUE),
    total_G = sum(G, na.rm = TRUE),
    total_V = sum(V, na.rm = TRUE),
    total_B = sum(B, na.rm = TRUE),
    mean_FAR = mean(FAR, na.rm = TRUE),
    median_FAR = median(FAR, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, output_csv, row.names = FALSE)

message("Saved explicit-allocation cells to: ", output_gpkg)
message("Saved summary table to: ", output_csv)
