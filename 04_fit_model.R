library(dplyr)
library(purrr)
library(broom)

 df_mod <- far_std_sf %>%
  filter(FAR > 0,
         distCBD > 0) %>%
  mutate(
    lnFAR = log(FAR),
    lnDist = log(distCBD)
  )


fit_city_models <- function(city_df) {
  
  # Linear: FAR ~ dist
  mod_linear <- tryCatch(
    lm(FAR ~ distCBD, data = city_df),
    error = function(e) NULL
  )
  
  # Exponential: log(FAR) ~ dist
  mod_exp <- tryCatch(
    lm(lnFAR ~ distCBD, data = city_df),
    error = function(e) NULL
  )
  
  # Gravitational: log(FAR) ~ log(dist)
  mod_grav <- tryCatch(
    lm(lnFAR ~ lnDist, data = city_df),
    error = function(e) NULL
  )
  
  # Log–log: log(FAR) ~ log(distCBD)
  mod_loglog <- tryCatch(
    lm(lnFAR ~ lnDist, data = city_df),
    error = function(e) NULL
  )
  
  tibble(
    model = c("linear", "exponential", "gravitational", "loglog"),
    fit = list(mod_linear, mod_exp, mod_grav, mod_loglog)
  )
}


results <- df_mod %>%
  group_by(city) %>%
  group_modify(~ fit_city_models(.x))

model_stats <- results %>%
  mutate(tidy = map(fit, tidy),
         glance = map(fit, glance)) %>%
  select(city, model, tidy, glance)


coef_table <- model_stats %>%
  unnest(tidy)

fit_stats <- model_stats %>%
  unnest(glance)

fit_and_summarise <- function(city_df, city_name) {
  
  models <- list(
    linear        = lm(FAR ~ distCBD, data = city_df),
    exponential   = lm(lnFAR ~ distCBD, data = city_df),
    gravitational = lm(lnFAR ~ lnDist,  data = city_df),
    loglog        = lm(lnFAR ~ lnDist,  data = city_df)
  )
  
  map_dfr(names(models), function(m) {
    mod <- models[[m]]
    
    # extract coef
    coefs <- tidy(mod)
    
    intercept <- coefs$estimate[coefs$term == "(Intercept)"]
    slope     <- coefs$estimate[coefs$term != "(Intercept)"][1]
    
    # extract model statistics
    gs <- glance(mod)
    
    tibble(
      City      = city_name,
      Intercept = intercept,
      Slope     = slope,
      R2        = gs$r.squared,
      n         = gs$df.residual + gs$df.null - gs$df.residual,
      Model     = m
    )
  })
}


results_df <- df_mod %>%
  group_by(city) %>%
  group_modify(~ fit_and_summarise(.x, .y$city)) %>%
  ungroup()

all_city_grids2 <- all_city_grids2 %>%
  mutate(
    dist_km = distCBD/1000,
    logFAR = log(FAR),
    logDist = log(dist_km)
  )
