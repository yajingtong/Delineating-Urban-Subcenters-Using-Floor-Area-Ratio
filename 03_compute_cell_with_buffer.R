library(sf)
library(dplyr)
library(purrr)

centers <- readRDS("data/centers.rds")
b_city <- readRDS("outputs/all_cities_buildings.rds")

process_city_buffers <- function(city_name,
                                 center_geom,
                                 buildings_sf,
                                 buffer_dist) {
  
  message("Processing city: ", city_name)
  
  #  Create grid
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
  
  #  Distance to CBD 
  cent <- st_coordinates(st_centroid(grid_sf))
  center_xy <- st_coordinates(center_geom)
  
  grid_sf <- grid_sf %>%
    mutate(
      X = cent[,1],
      Y = cent[,2],
      distCBD = sqrt((X - center_xy[1])^2 + (Y - center_xy[2])^2),
      city = city_name
    )
  
  #  Select city buildings
  b_city <- buildings_sf %>% filter(city == city_name)
  
  if (nrow(b_city) == 0) {
    grid_sf$buffer_area <- NA_real_
    return(grid_sf)
  }
  
  #   Buffer buildings (NO union)  
  b_buffered <- st_buffer(b_city, dist = buffer_dist)
  
  #   Crop buffers to grid cells 
  buffers_in_grids <- st_intersection(
    grid_sf %>% select(id),
    b_buffered
  )
  
  #   Aggregate buffer area per grid 
  buffer_stats <- buffers_in_grids %>%
    mutate(buf_area = st_area(geometry)) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(
      buffer_area = sum(as.numeric(buf_area), na.rm = TRUE),
      n_buffers   = n(),
      .groups = "drop"
    )
  
  #   Join back to grid  
  grid_sf <- grid_sf %>%
    left_join(buffer_stats, by = "id")
  
  return(grid_sf)
}



all_city_grids <- map_dfr(
  1:nrow(centers),
  function(i) {
    process_city_buffers(
      city_name   = centers$Name[i],
      center_geom = centers[i, ],
      buildings_sf = b_city,      
      buffer_dist = buffer_dist
    )
  }
)




process_city_grid <- function(city_name,
                              all_city_grids,
                              buildings_sf,
                              buffer_dist) {
  
  message("Processing city: ", city_name)
  
  # grid
  grid_city <- all_city_grids %>% filter(city == city_name)
  
  # buildings 
  b_city <- buildings_sf %>% filter(city == city_name)
  
  # No buildings
  if (nrow(b_city) == 0) {
    grid_city <- grid_city %>%
      mutate(
        buffer_area = NA_real_,
        G = NA_real_,
        V = NA_real_,
        H = NA_real_,
        N = 0L,
        FAR = NA_real_
      )
    
    return(grid_city)
  }
  
  
  # buffer without union
  b_buffered <- b_city %>%
    st_buffer(dist = buffer_dist)
  
  # intersect with grid
  inter <- st_intersection(
    grid_city %>% select(id),
    b_buffered
  )
  
  if (nrow(inter) == 0) {
    grid_city <- grid_city %>%
      mutate(
        buffer_area = 0,
        G = 0,
        V = 0,
        H = NA_real_,
        N = 0L,
        FAR = NA_real_
      )
    
    
    return(grid_city)
  }
  

  
  # aggregate per grid
  stats <- inter %>%
    mutate(buf_area = st_area(geometry)) %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(
      buffer_area = sum(as.numeric(buf_area), na.rm = TRUE),
      G = sum(gi, na.rm = TRUE),
      V = sum(vi, na.rm = TRUE),
      H = mean(hi, na.rm = TRUE),
      B = sum(bi,na.rm = TRUE),
      N = n(), # number of buffers intersecting
      .groups = "drop"
    ) %>%
    mutate(
      FAR = if_else(buffer_area > 0, G / buffer_area, NA_real_)
    )
  
 
  # ---- join back ----
  grid_city <- grid_city %>%
    left_join(stats, by = "id")
  
  return(grid_city)
}

cities <- unique(all_city_grids_buffered$city)

results_list <- vector("list", length(cities))
names(results_list) <- cities

for (i in seq_along(cities)) {
  
  ct <- cities[i]
  
  res <- process_city_grid(
    city_name      = ct,
    all_city_grids = all_city_grids,
    buildings_sf   = b_city,
    buffer_dist    = buffer_dist
  )
  
  results_list[[i]] <- res
}


all_city_2 <- do.call(rbind, results_list)



plot(all_city_grids_buffered2$distCBD ~ all_city_grids_buffered2$FAR)








process_city_grid <- function(city_name,
                              all_city_grids,
                              buildings_sf,
                              buffer_dist) {
  
  message("Processing city: ", city_name)
  
  # grid
  grid_city <- all_city_grids %>% 
    filter(city == city_name)
  
  # buildings
  b_city <- buildings_sf %>% 
    filter(city == city_name)
  
  # no buildings
  if (nrow(b_city) == 0) {
    return(
      grid_city %>%
        mutate(
          buffer_area = NA_real_,
          G = NA_real_,
          V = NA_real_,
          H = NA_real_,
          N = 0L,
          FAR = NA_real_
        )
    )
  }
  
  # buffer buildings (NO union here)
  b_buffered <- b_city %>%
    st_buffer(dist = buffer_dist)
  
  # intersect grid × buffers
  inter <- st_intersection(
    grid_city %>% select(id),
    b_buffered
  )
  
  if (nrow(inter) == 0) {
    return(
      grid_city %>%
        mutate(
          buffer_area = 0,
          G = 0,
          V = 0,
          H = NA_real_,
          N = 0L,
          FAR = NA_real_
        )
    )
  }
  
  # union buffers PER CELL (geometry only) 
  buffer_area_df <- inter %>%
    group_by(id) %>%
    summarise(
      geometry = st_union(geometry),
      .groups = "drop"
    ) %>%
    mutate(
      buffer_area = as.numeric(st_area(geometry))
    ) %>%
    st_drop_geometry()
  
  # ---- 2. aggregate attributes PER CELL ----
  attr_df <- inter %>%
    st_drop_geometry() %>%
    group_by(id) %>%
    summarise(
      G = sum(gi, na.rm = TRUE),
      V = sum(vi, na.rm = TRUE),
      H = mean(hi, na.rm = TRUE),
      N = n(),
      .groups = "drop"
    )
  
  stats <- buffer_area_df %>%
    left_join(attr_df, by = "id") %>%
    mutate(
      FAR = if_else(buffer_area > 0, G / buffer_area, NA_real_)
    )
  
  grid_city <- grid_city %>%
    left_join(stats, by = "id")
  
  return(grid_city)
}



cities <- unique(all_city_grids$city)

results_list <- vector("list", length(cities))
names(results_list) <- cities

for (i in seq_along(cities)) {
  
  ct <- cities[i]
  
  results_list[[i]] <- process_city_grid(
    city_name      = ct,
    all_city_grids = all_city_grids_buffered,
    buildings_sf   = b_city,
    buffer_dist    = buffer_dist
  )
}

all_city_grids2<- do.call(rbind, results_list)






library(ggplot2)





ggplot(
  all_city_grids2,
  aes(x = distCBD, y = FAR)
) +
  geom_point(alpha = 0.05) +
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = TRUE
  ) +
  labs(
    x = "Distance to CBD",
    y = "FAR",
    title = "Global FAR–Distance Relationship (All Cities Pooled)"
  ) +
  theme_minimal()



library(dplyr)


far_km_city <- all_city_grids2%>%
  filter(!is.na(FAR), distCBD <=20000) %>%
  mutate(
    dist_km = floor(distCBD / 1000)  # use floor(distCBD) if already in km
  ) %>%
  group_by(city, dist_km) %>%
  summarise(
    FAR_mean = mean(FAR, na.rm = TRUE),
    n_cells  = n(),
    .groups = "drop"
  ) # optional but recommended


label_df <- far_km_city %>%
  group_by(city) %>%
  filter(dist_km == max(dist_km, na.rm = TRUE)) %>%
  ungroup()



library(plotly)

p <- plot_ly()

# one trace per city
for (ct in unique(far_km_city$city)) {
  
  df_ct <- far_km_city %>% filter(city == ct)
  
  p <- p %>%
    add_lines(
      data = df_ct,
      x = ~dist_km,
      y = ~FAR_mean,
      name = ct,
      hovertemplate =
        paste0(
          "<b>", ct, "</b><br>",
          "Distance: %{x} km<br>",
          "log(Mean FAR): %{y:.3f}<extra></extra>"
        ),
      opacity = 0.5
    )
}

# add city labels at curve ends
p <- p %>%
  add_text(
    data = label_df,
    x = ~dist_km,
    y = ~FAR_mean,
    text = ~city,
    textposition = "middle right",
    showlegend = FALSE
  )

p <- p %>%
  layout(
    title = "Mean FAR per km as a Function of Distance to CBD",
    xaxis = list(title = "Distance to CBD (km)"),
    yaxis = list(title = "Mean FAR"),
    hovermode = "closest"
  )

p

ggplot(
  data = all_city_grids2,
  mapping = aes(x = dist_km, y = FAR_mean, color = island)
) +
  geom_point() +
  geom_smooth(se = FALSE)




p_decay <- ggplot(plot_df, aes(x = logDist, y = logFAR)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = lm, color = "red", fill = "#69b3a2", se = TRUE) +
  facet_wrap(~ city, scales = "free") +
  labs(
    x = "log(Distance to CBD) [km]",
    y = "log(FAR)",
    title = "Distance Decay of FAR per City",
    subtitle = "Linear fit in log–log space (gravitational form)"
  ) +
  theme()

p_decay



p_decay2 <- ggplot(plot_df, aes(x = distCBD/1000, y = logFAR)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = lm, color = "red", fill = "#69b3a2", se = TRUE) +
  facet_wrap(~ city, scales = "free_x") +
  labs(
    x = "Distance to CBD (km)",
    y = "log(FAR)"
  ) +
  theme()

p_decay2



df_model <- all_city_grids2%>%
  filter(FAR > 0, distCBD > 0) %>%
  mutate(dist_km = distCBD / 1000)

fit_grav_city <- function(df) {
  
  # get starting values from log–log linear fit
  start_lm <- lm(log(FAR) ~ log(dist_km), data = df)
  a_start  <- exp(coef(start_lm)[1])
  b_start  <- coef(start_lm)[2]
  
  nls(
    FAR ~ a * dist_km^b,
    data = df,
    start = list(a = a_start, b = b_start),
    control = nls.control(maxiter = 200, warnOnly = TRUE)
  )
}


grav_models <- df_model %>%
  group_split(city) %>%
  setNames(unique(df_model$city)) %>%
  map(fit_grav_city)

pred_df <- map2_df(
  grav_models,
  names(grav_models),
  function(mod, city_name) {
    
    d_seq <- seq(
      min(df_model$dist_km[df_model$city == city_name]),
      max(df_model$dist_km[df_model$city == city_name]),
      length.out = 200
    )
    
    tibble(
      city = city_name,
      dist_km = d_seq,
      FAR_pred = predict(mod, newdata = list(dist_km = d_seq))
    )
  }
)


p_grav <- ggplot(df_model, aes(x = dist_km, y = FAR)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_line(data = pred_df, aes(x = dist_km, y = FAR_pred),
            color = "red", linewidth = 1) +
  facet_wrap(~ city, scales = "free_y") +
  labs(
    x = "Distance to CBD (km)",
    y = "FAR",
    title = "Nonlinear Gravitational Decay of FAR",
    subtitle = "FAR = α · Distance^β"
  ) +
  theme_minimal()



p_grav



library(dplyr)
library(sf)
library(ggplot2)
library(viridis)

far_std_sf <- all_city_grids2 %>%  filter(dist_km <= 5) 

far_std_sf <- far_std_sf %>%  
  group_by(city) %>% 
  mutate(
    FAR_std = (FAR - min(FAR, na.rm = TRUE)) /
      (max(FAR, na.rm = TRUE) - min(FAR, na.rm = TRUE))
  ) %>%
  ungroup()




p_map <- ggplot(far_std_sf) +
  geom_sf(aes(fill = FAR_std), color = NA) +
  scale_fill_viridis(
    option = "viridis",
    name = "Standardized FAR",
    limits = c(0, 1)
  ) +
  facet_wrap(~ city) +  
  labs(
    title = "Standardized FAR Distribution Across Cities",
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "right"
  )
library(purrr)

city_fits <- far_std_sf %>%
  group_by(city) %>%
  group_modify(~{
    
    m <- lm(log(FAR) ~ log(dist_km), data = .x)
    
    tibble(
      alpha = exp(coef(m)[1]),
      beta  = coef(m)[2],
      R2    = summary(m)$r.squared
    )
  }) %>%
  ungroup()

curve_df <- city_fits %>%
  group_by(city) %>%
  group_modify(~{
    
    d_seq <- seq(
      min(df_model$dist_km[df_model$city == .y$city]),
      max(df_model$dist_km[df_model$city == .y$city]),
      length.out = 200
    )
    
    tibble(
      dist_km = d_seq,
      FAR_fit = .x$alpha * d_seq^.x$beta
    )
  }) %>%
  ungroup()



ggplot(df_model, aes(x = dist_km, y = FAR)) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_line(
    data = curve_df,
    aes(y = FAR_fit),
    color = "red",
    linewidth = 1
  ) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ city, scales = "free_y") +
  labs(
    title = "City-Specific Gravitational Fits of FAR",
    subtitle = "FAR = α · distance^β",
    x = "Distance to CBD (km)",
    y = "FAR"
  ) +
  theme_minimal()




library(purrr)

city_models_expgrav <- far_std_sf %>%
  group_split(city) %>%
  setNames(unique(far_std_sf$city)) %>%
  map(fit_exp_grav)

fit_exp_grav <- function(df) {
  
  # starting values from log-log gravity
  lm_start <- lm(log(FAR) ~ log(dist_km), data = df)
  
  a0 <- exp(coef(lm_start)[1])
  b0 <- coef(lm_start)[2]
  g0 <- 0.05   # weak exponential cutoff
  
  nls(
    FAR ~ a * dist_km^b * exp(-g * dist_km),
    data = df,
    start = list(a = a0, b = b0, g = g0),
    algorithm = "port",
    lower = c(a = 0, b = -5, g = 0),
    control = nls.control(maxiter = 500, warnOnly = TRUE)
  )
}


curve_df <- map2_df(
  city_models_expgrav,
  names(city_models_expgrav),
  function(mod, city_name) {
    
    df_sub <- df_city %>% filter(city == city_name)
    
    d_seq <- seq(
      min(df_sub$dist_km),
      max(df_sub$dist_km),
      length.out = 200
    )
    
    tibble(
      city = city_name,
      dist_km = d_seq,
      FAR_fit = predict(mod, newdata = list(dist_km = d_seq))
    )
  }
)


library(ggplot2)

ggplot(df_city, aes(x = dist_km, y = FAR)) +
  geom_point(alpha = 0.3, size = 0.7) +
  geom_line(
    data = curve_df,
    aes(y = FAR_fit),
    color = "red",
    linewidth = 1
  ) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~ city, scales = "free_y") +
  labs(
    title = "Exponential–Gravity Distance Decay of FAR",
    subtitle = expression(FAR == alpha %.% d^beta %.% e^{-gamma %.% d}),
    x = "Distance to CBD (km)",
    y = "FAR"
  ) +
  theme_minimal()

city_models <- far_std_sf %>%
  group_by(city) %>%
  group_map(~ lm(log(FAR) ~ log(dist_km), data = .x)) %>%
  setNames(unique(df_city$city))



city_params <- map_dfr(
  city_models,
  ~ tibble(
    alpha = exp(coef(.x)[1]),
    beta  = coef(.x)[2]
  ),
  .id = "city"
)


curve_df <- far_std_sf %>%
  left_join(df_city %>% group_by(city) %>%
              summarise(
                dmin = min(dist_km),
                dmax = max(dist_km)
              ),
            by = "city") %>%
  rowwise() %>%
  mutate(
    data = list(
      tibble(
        dist_km = seq(dmin, dmax, length.out = 200),
        FAR_fit = alpha * dist_km^beta
      )
    )
  ) %>%
  unnest(data)





# 6. Plotting
ggplot() +
  # Individual City Lines (Power fits)
  geom_line(data = plot_data, 
            aes(x = DIST, y = FAR, group = city.x, color = Pop2025), 
            size = 1, alpha = 0.6) +
  # General Model Line (Highlighted)
  geom_line(data = general_curve, 
            aes(x = DIST, y = FAR), 
            color = "black", size = 1.5, linetype = "dashed") +
  # Log-Log Scales (Makes the curves appear as straight lines)
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50)) +
  scale_y_log10(labels = comma) +
  scale_color_viridis_c(option = "magma", direction = -1, labels = comma) +
  # Labels
  labs(
    title = "Log-Log Gravitational Model: FAR vs DIST",
    subtitle = "Colored by Population (2025); Black Dashed Line = General Model",
    x = "Log Distance (km)",
    y = "Log Floor Area Ratio (FAR)",
    color = "Population"
  ) +
  theme_minimal() +
  theme(legend.position = "right")













pl_fits <- far_std_sf %>%
  filter(FAR > 0, dist_km > 0) %>%
  group_by(city) %>%
  do(
    tidy(
      lm(log(FAR) ~ log(dist_km), data = .)
    )
  ) %>%
  ungroup()



pl_params <- pl_fits %>%
  select(city, term, estimate) %>%
  tidyr::pivot_wider(
    names_from = term,
    values_from = estimate
  ) %>%
  rename(
    log_a = `(Intercept)`,
    b = `log(dist_km)`
  ) %>%
  mutate(a = exp(log_a))


pred_data <- far_std_sf %>%
  group_by(city) %>%
  summarise(
    dist_min = min(dist_km, na.rm = TRUE),
    dist_max = max(dist_km, na.rm = TRUE)
  ) %>%
  left_join(pl_params, by = "city") %>%
  rowwise() %>%
  do({
    tibble(
      city = .$city,
      dist_km = seq(.$dist_min, .$dist_max, length.out = 100),
      FAR_fit = .$a * dist_km ^ .$b
    )
  }) %>%
  ungroup()

dist_seq <- seq(
  min(far_std_sf$dist_km, na.rm = TRUE),
  max(far_std_sf$dist_km, na.rm = TRUE),
  length.out = 300
)

general_fits <- tibble(
  dist_km = dist_seq,
  FAR_linear = 1.0199009 + (-4.676987e-05) * dist_seq,
  FAR_exp    = 0.6303637 * exp(-8.849348e-05 * dist_seq),
  FAR_power  = 6.6892177 * dist_seq ^ (-9.043024e-01)
)

ggplot() +
  geom_line(
    data = pred_data,
    aes(dist_km, FAR_fit, color = city),
    linewidth = 0.8
  ) +
  
  # general linear
  geom_line(
    data = general_fits,
    aes(dist_km, FAR_linear),
    linewidth = 1.4,
    linetype = "dashed",
    color = "black"
  ) +
  
  # general negative exponential
  geom_line(
    data = general_fits,
    aes(dist_km, FAR_exp),
    linewidth = 1.4,
    linetype = "dotdash",
    color = "black"
  ) +
  
  geom_line(
    data = general_fits,
    aes(dist_km, FAR_power),
    linewidth = 1.8,
    linetype = "solid",
    color = "black"
  ) + scale_color_viridis_d(option = "viridis", direction = -1) + 

  labs(
    x = "Distance to CBD (km)",
    y = "FAR",
    #title = "City-specific power-law fits with global decay models",
    subtitle = "Black lines: Linear (dashed), Negative exponential (dotdash), Power law (solid)"
  ) +
  
  theme_minimal()



ggplot(plot_df, aes(x = dist_km, y = FAR, color = city)) +
  scale_color_viridis_d(option = "viridis", direction = -1) + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(title = "Standardised Floor area ratio gradients by city deceding with population",
       x = "Distance to CBD (km)",
       y = "Standardised FAR") 
theme_minimal()




ggplot() +
  geom_line(
    data = pred_data,
    aes(x = dist_km, y = FAR_fit, color = city),
    linewidth = 1
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Distance to CBD (km)",
    y = "FAR",
    #title = "Log-log fits "
  ) +
  theme_minimal()





ggplot() +
  geom_line(
    data = pred_data,
    aes(dist_km, FAR_fit, color = city),
    linewidth = 0.8
  ) + scale_color_viridis_d(option = "viridis", direction = -1) + 
  labs(
    x = "Distance to CBD (km)",
    y = "FAR",
    #title = "City-specific power-law fits with global decay models",
  ) +
  theme_minimal()


pred_data <- pred_data  %>%
  mutate(city = reorder(city, -Pop2025)) 


pred_data <- pred_data %>%
     left_join(plot_data, by = "city")


