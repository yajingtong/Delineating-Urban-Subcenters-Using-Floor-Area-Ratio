
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(scales)

input_csv <- paste0(
  "outputs/far_power_law_positive_residual_hh_patches_41/",
  "urban_configuration_typology_3km/",
  "city_urban_configuration_centroid_distance_3km_41_cities.csv"
)
output_dir <- "outputs/zipf_law_all_41_cities"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

city_data <- read_csv(input_csv, show_col_types = FALSE) |>
  transmute(
    city = as.character(city),
    population = as.numeric(population_2025)
  ) |>
  filter(!is.na(city), is.finite(population), population > 0) |>
  distinct(city, .keep_all = TRUE) |>
  arrange(desc(population), city) |>
  mutate(
    rank = row_number(),
    zipf_population = first(population) / rank
  )

if (nrow(city_data) != 41) {
  warning("Expected 41 cities but found ", nrow(city_data), ".")
}

rank_size_model <- lm(log10(population) ~ log10(rank), data = city_data)
model_summary <- summary(rank_size_model)
intercept <- unname(coef(rank_size_model)[1])
slope <- unname(coef(rank_size_model)[2])

city_data <- city_data |>
  mutate(fitted_population = 10^(intercept + slope * log10(rank)))

fit_stats <- tibble(
  cities = nrow(city_data),
  intercept_log10 = intercept,
  slope = slope,
  slope_p_value = coef(model_summary)["log10(rank)", "Pr(>|t|)"],
  r_squared = model_summary$r.squared,
  adjusted_r_squared = model_summary$adj.r.squared,
  theoretical_zipf_slope = -1
)

plot_subtitle <- sprintf(
  "Fitted slope = %.3f; R² = %.3f | theoretical Zipf slope = -1",
  slope,
  model_summary$r.squared
)

p <- ggplot(city_data, aes(x = rank, y = population)) +
  geom_line(
    aes(y = zipf_population, linetype = "Theoretical Zipf law"),
    color = "#C7432B",
    linewidth = 0.85
  ) +
  geom_line(
    aes(y = fitted_population, linetype = "Fitted rank-size relationship"),
    color = "#163D5C",
    linewidth = 0.9
  ) +
  geom_point(
    aes(fill = log10(population)),
    shape = 21,
    size = 3.0,
    stroke = 0.55,
    color = "#142D3D"
  ) +
  geom_text_repel(
    aes(label = city),
    size = 2.35,
    color = "#273B48",
    seed = 2026,
    box.padding = 0.28,
    point.padding = 0.22,
    min.segment.length = 0,
    segment.color = "#8BA2AF",
    segment.size = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_x_log10(
    breaks = c(1, 2, 3, 5, 10, 20, 30, 41),
    labels = label_number(accuracy = 1),
    expand = expansion(mult = c(0.035, 0.14))
  ) +
  scale_y_log10(
    labels = label_number(scale = 1e-6, suffix = " M", accuracy = 0.1),
    expand = expansion(mult = c(0.24, 0.18))
  ) +
  scale_fill_gradientn(
    colours = c("#DCEEF4", "#69A9C4", "#163D5C"),
    guide = "none"
  ) +
  scale_linetype_manual(
    values = c(
      "Fitted rank-size relationship" = "solid",
      "Theoretical Zipf law" = "22"
    ),
    breaks = c("Theoretical Zipf law", "Fitted rank-size relationship"),
    name = NULL
  ) +
  labs(
    title = "Zipf's law across 41 French urban areas",
    subtitle = plot_subtitle,
    x = "Population rank (log scale)",
    y = "Population in 2025 (log scale)",
    caption = "The theoretical line is anchored to the population of the largest urban area."
  ) +
  theme_minimal(base_family = "sans", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = "#173A4F"),
    plot.subtitle = element_text(size = 9.5, color = "#4B6170"),
    plot.caption = element_text(size = 8, color = "#657681", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#D9E3E8", linewidth = 0.35),
    axis.title = element_text(face = "bold", color = "#274553"),
    legend.position = "inside",
    legend.position.inside = c(0.72, 0.86),
    legend.background = element_rect(fill = alpha("white", 0.9), color = "#C7D6DE"),
    legend.key.width = grid::unit(1.5, "cm"),
    plot.margin = margin(10, 26, 14, 10)
  )

write_csv(
  city_data,
  file.path(output_dir, "zipf_law_all_41_cities_plot_data.csv")
)
write_csv(
  fit_stats,
  file.path(output_dir, "zipf_law_all_41_cities_fit_statistics.csv")
)

ggsave(
  file.path(output_dir, "zipf_law_all_41_cities.png"),
  p,
  width = 18,
  height = 15,
  units = "cm",
  dpi = 400,
  bg = "white"
)
ggsave(
  file.path(output_dir, "zipf_law_all_41_cities.pdf"),
  p,
  width = 18,
  height = 15,
  units = "cm",
  device = pdf,
  useDingbats = FALSE,
  bg = "white"
)

message("Saved Zipf plot and model results to: ", output_dir)
