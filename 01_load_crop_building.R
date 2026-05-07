library(sf)
library(purrr)
library(dplyr)

#For each city centers:
#Compute its location.
#Loop over all building shapefiles, read their bounding box centroid.
#Allocate the building file is closest to that city centers.
#Load the shapefile, assign each file with city name, and store it.

# Paths
# Read city centers
centers <- readRDS("data/centers.rds")
target_crs <- 32632
centers <- st_transform(centers, target_crs)

# Parameters
width_ring <- 250 
max_rings <- 100
radius_max <- width_ring * max_rings # 25 km radius

# List all building tiles
building_files <- list.files("data/building/fr",pattern = "\\.shp$",recursive = TRUE,full.names = TRUE)

# Precompute bbox for each tile 
tile_info <- map(building_files, function(f) {
  tmp <- st_read(f, quiet = TRUE)   # only reads header + bbox
  bb  <- st_as_sfc(st_bbox(tmp), crs = st_crs(tmp))
  bb  <- st_transform(bb, target_crs)
  list(file = f, bbox = bb)})


all_cities_buildings <- map_dfr(1:nrow(centers), function(i) {
  
  city_name <- centers$Name[i]
  city_geom <- centers[i, ]
  
  message("\n=== Processing: ", city_name, " ===")
  
  #create  buffer per city
  city_buffer <- st_buffer(city_geom, radius_max)
  
  # Find all building tiles overlapping the buffer
  
  overlaps <- sapply(tile_info, function(ti)
    lengths(st_intersects(ti$bbox, city_buffer)) > 0
  )
  
  tile_subset <- tile_info[overlaps]
  
  message("Found tiles: ", length(tile_subset))
  
  # Stop if none tile found
  if (length(tile_subset) == 0) {
    warning("No building tiles found for ", city_name)
    return(NULL)
  }
  
   #  Load, crop, and combine all matching tiles
  b_city <- map_dfr(tile_subset, function(ti) {
    
    f <- ti$file
    message("  → Loading tile: ", basename(f))
    
    b <- st_read(f, quiet = TRUE, options = "PROMOTE_TO_MULTI=YES")
    
    # match CRS
    if (st_crs(b) != st_crs(target_crs)) {
      b <- st_transform(b, target_crs)
    }
    
    # Check&Fix invalid geometries  
    if (any(!st_is_valid(b))) {
      b <- st_make_valid(b)
    }
    
    # Crop using intersects 
    idx <- st_intersects(b, city_buffer) %>% lengths() > 0
    b_crop <- b[idx, ]
    b_crop
  })
  
  b_city$city <- city_name
  
  message("  Total buildings kept for ", city_name, ": ", nrow(b_city))
  
  b_city #return building sf
})

#compute building metrics
all_cities_buildings <- all_cities_buildings %>%
  mutate(
    bi = as.numeric(st_area(geometry)),       # footprint (m²)
    hi = as.numeric(Height),                  # height (m)
    si = pmax(1, round(hi / 3)),              # storeys
    gi = bi * si,                             # floor area (m²)
    vi = bi * hi ) 

#save building sf
st_write(all_cities_buildings,"outputs/all_city.gpkg")
building_points <- st_centroid(all_cities_buildings)
saveRDS(building_points,"outputs/buildings_points.rds")
