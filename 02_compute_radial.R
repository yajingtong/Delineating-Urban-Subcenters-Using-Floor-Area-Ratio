library(sf)
library(dplyr)
library(purrr)
library(future)
library(lubridate)

plan(multisession, workers = 8)

ComputeRadialProfiles_Buildings <- function(buildings, centers,
                                            no_column_name = "Name",
                                            width_ring = 250,
                                            max_rings = 100,
                                            floor_height = 3,
                                            buffer_list = c(50, 100, 200),
                                            output_dir = "RadialProfiles_Results") {
  
  all_start <- Sys.time()
  city_names <- unique(buildings$city)
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  for (city in city_names) {
    city_start <- Sys.time()
    message("Processing: ", city)

    #  city data
    city_center <- centers %>% filter(!!sym(no_column_name) == city)
    b_city <- buildings %>% filter(city == !!city)
    
    if (nrow(b_city) == 0 | nrow(city_center) == 0) {
      next
    }
    
    city_center <- st_transform(city_center, st_crs(b_city))
    
    #   Clip to max analysis extent 
    city_boundary <- st_buffer(city_center, dist = max_rings * width_ring)
    b_city <- st_intersection(b_city, city_boundary)
    
    #  Compute distances and metrics per building 
    centroids <- st_centroid(b_city)
    distances <- as.numeric(st_distance(centroids, city_center))
    
    b_city <- b_city %>%
      mutate(
        ri = ceiling(distances / width_ring),     # ring number
        bi = as.numeric(st_area(geometry)),       # footprint (m²)
        hi = as.numeric(Height),                  # height (m)
        si = pmax(1, round(hi / floor_height)),   # storeys
        gi = bi * si,                             # floor area (m²)
        vi = bi * hi                              # volume (m³)
      )
    
    #  Build ring polygons 
    rings_sf <- map_dfr(1:max_rings, function(i) {
      outer <- st_buffer(city_center, dist = i * width_ring)
      inner <- if (i > 1) st_buffer(city_center, dist = (i - 1) * width_ring) else NULL
      ring_geom <- if (!is.null(inner)) st_difference(outer, inner) else outer
      st_sf(ri = i, geometry = st_geometry(ring_geom))
    })
    rings_sf <- st_make_valid(rings_sf)
    rings_sf$Ar <- as.numeric(st_area(rings_sf))
    
    #   Aggregate building metrics per ring (without buffer)  
    base_data <- b_city %>%
      st_drop_geometry() %>%
      group_by(ri) %>%
      summarise(
        xr = n(),                          # number of buildings
        Br = sum(bi, na.rm = TRUE),        # total footprint
        Gr = sum(gi, na.rm = TRUE),        # total floor area
        Hr = sum(hi, na.rm = TRUE),        # total height
        Vr = sum(vi, na.rm = TRUE),        # total volume
        br = mean(bi, na.rm = TRUE),       # average footprint
        gr = mean(gi, na.rm = TRUE),       # average floor area 
        hr = mean(hi, na.rm = TRUE),       # average height
        .groups = "drop"
      ) %>%
      left_join(rings_sf %>% st_drop_geometry() %>% select(ri, Ar), by = "ri") %>%
      mutate(
        Xr  = xr / (Ar / 1e6),
        Br  = Br / 1e6,
        city = city
      )
    
    #  Loop through multiple buffer distances  
    for (buffer_dist in buffer_list) {
      message("   Computing buffer coverage for buffer = ", buffer_dist, " m")
      
      # Create dissolved buffers
      b_buffered <- b_city %>%
        st_buffer(dist = buffer_dist) %>%
        st_union() %>%
        st_sf(geometry = .)
      
      # Intersect buffer with rings
      buffers_in_rings <- st_intersection(rings_sf, b_buffered)
      buffers_in_rings$A_buffer <- as.numeric(st_area(buffers_in_rings))
      
      # Compute ratios per ring
      ring_ratios <- buffers_in_rings %>%
        st_drop_geometry() %>%
        group_by(ri) %>%
        summarise(A_buffer = sum(A_buffer, na.rm = TRUE), .groups = "drop") %>%
        left_join(rings_sf %>% st_drop_geometry() %>% select(ri, Ar), by = "ri") %>%
        mutate(Rr = A_buffer / Ar)
      
      # Total buffer/ring area ratio
      R_total <- sum(ring_ratios$A_buffer, na.rm = TRUE) / sum(ring_ratios$Ar, na.rm = TRUE)
      
      # Combine with  radial metrics
      Synth_data <- base_data %>%
        left_join(ring_ratios %>% select(ri, Rr,A_buffer), by = "ri") %>%
        mutate(
          R_total = R_total,
          buffer_m = buffer_dist,
          FAR = Gr / A_buffer
        )
      
      # Save each buffer result
      out_file <- file.path(output_dir, paste0("RadialProfiles_", city, "_buf", buffer_dist, "m.csv"))
      write.csv(Synth_data, out_file, row.names = FALSE)}
    
    
    city_end <- Sys.time()
    message(" finished ", city, " in ", round(difftime(city_end, city_start, units = "mins"), 2), " min")
  }
  
  all_end <- Sys.time()
  message("\nAll cities processed in ", 
          round(difftime(all_end, all_start, units = "mins"), 2), " minutes ")}



# call function ComputeRadialProfiles_Buildings()
ComputeRadialProfiles_Buildings(
  buildings = all_cities_buildings,
  centers   = centers,
  no_column_name = "Name",
  width_ring = 250,
  max_rings  = 100,
  floor_height = 3,
  buffer_list = c(50),
  output_dir = "Radial-Profiles")
  #buffer_list = c(50, 100, 150,200,250,300,350,400,450,500,550,600))



