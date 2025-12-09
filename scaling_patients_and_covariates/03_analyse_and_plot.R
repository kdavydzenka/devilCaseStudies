
rm(list = ls())
library(tidyverse)
library(ggh4x)
library(patchwork)
library(scales)

# Define color palette
MY_PALETTE <- c(
  "Devil (base)" = "#099668", "Devil (mixed)" = "#099668", 
  "Devil" = "#099668", "devil" = "#099668",
  "Nebula" = "steelblue", "NEBULA" = "steelblue",
  "edgeR" = "#7D629E", "edgeR (Pb)" = "#7D629E",
  "limma" = "#B96461", "limma (Pb)" = "#B96461",
  "glmGamPoi (cell)" = "#EAB578", "glmGamPoi" = "#EAB578",
  "limmaDupCorr (cell)" = "#8B0000", "limmaDupCorr" = "#8B0000",
  "Seurat (cell)" = "#708090", "Seurat (Wilcox)" = "#708090", "Seurat" = "#708090",
  "MAST (cell)" = "#D8BFD8", "MAST" = "#D8BFD8"
)

# Load and prepare data
res <- readRDS("summarized_results/results.rds") %>%
  na.omit() %>%
  pivot_longer(
    cols = !c(ncells, n_patients, n_covariates, n_genes, seed, continuous, idx),
    names_to = "method",
    values_to = "runtime"
  )

# Set base method for speedup calculations
BASE_METHOD <- "NEBULA"

# Figure 1: Scalability by Number of Patients ####

# Linear Model Fits (for global slope)
lm_fits <- res %>%
  group_by(method) %>%
  do({
    # Fit linear model on log-transformed runtime
    model <- lm(log10(runtime) ~ n_patients, data = .)
    data.frame(
      slope = coef(model)[2],
      intercept = coef(model)[1],
      r_squared = summary(model)$r.squared,
      p_value = summary(model)$coefficients[2, 4]
    )
  }) %>%
  ungroup() %>%
  arrange(desc(slope))

# LOESS-based slope measures
loess_avg_slope <- res %>%
  group_by(method) %>%
  do({
    # Fit LOESS model
    loess_model <- loess(log10(runtime) ~ n_patients, data = ., span = 0.75)
    
    # Predict across patient range
    patient_range <- seq(min(.$n_patients), max(.$n_patients), length.out = 100)
    predictions <- predict(loess_model, newdata = data.frame(n_patients = patient_range))
    
    # Calculate derivatives (differences)
    diffs <- diff(predictions) / diff(patient_range)
    
    data.frame(
      avg_slope = mean(diffs, na.rm = TRUE),
      median_slope = median(diffs, na.rm = TRUE),
      max_slope = max(diffs, na.rm = TRUE),
      min_slope = min(diffs, na.rm = TRUE)
    )
  }) %>%
  ungroup() %>%
  arrange(desc(avg_slope))

cat("\n=== LOESS Average Derivatives (local slopes) ===\n")
print(loess_avg_slope, n = Inf)

df_slopes_patients = dplyr::bind_rows(
  lm_fits %>% dplyr::select(method, slope) %>% dplyr::mutate(regression = "Linear"),
  loess_avg_slope %>% dplyr::select(method, avg_slope) %>% dplyr::rename(slope = avg_slope) %>% dplyr::mutate(regression = "LOESS")
) %>% dplyr::mutate(Variable = "Patients") %>% tidyr::pivot_wider(values_from = slope, names_from = regression, names_prefix = "Slope ") 

# Actual plot
p1_patients <- res %>%
  ggplot(aes(x = n_patients, y = runtime, color = method)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_y_continuous(transform = "log10") +
  scale_color_manual(values = MY_PALETTE, name = "Method") +
  labs(
    x = "Number of Patients",
    y = "Runtime (seconds, log scale)"
  ) +
  theme_bw()
p1_patients


# Figure 2: Scalability by Number of Covariates ####
lm_fits <- res %>%
  group_by(method) %>%
  do({
    # Fit linear model on log-transformed runtime
    model <- lm(log10(runtime) ~ n_covariates, data = .)
    data.frame(
      slope = coef(model)[2],
      intercept = coef(model)[1],
      r_squared = summary(model)$r.squared,
      p_value = summary(model)$coefficients[2, 4]
    )
  }) %>%
  ungroup() %>%
  arrange(desc(slope))

# LOESS-based slope measures

loess_avg_slope <- res %>%
  group_by(method) %>%
  do({
    # Fit LOESS model
    loess_model <- loess(log10(runtime) ~ n_covariates, data = ., span = 0.75)
    
    # Predict across patient range
    cov_range <- seq(min(.$n_covariates), max(.$n_covariates), length.out = 100)
    predictions <- predict(loess_model, newdata = data.frame(n_covariates = cov_range))
    
    # Calculate derivatives (differences)
    diffs <- diff(predictions) / diff(cov_range)
    
    data.frame(
      avg_slope = mean(diffs, na.rm = TRUE),
      median_slope = median(diffs, na.rm = TRUE),
      max_slope = max(diffs, na.rm = TRUE),
      min_slope = min(diffs, na.rm = TRUE)
    )
  }) %>%
  ungroup() %>%
  arrange(desc(avg_slope))

df_slopes_covariates = dplyr::bind_rows(
  lm_fits %>% dplyr::select(method, slope) %>% dplyr::mutate(regression = "Linear"),
  loess_avg_slope %>% dplyr::select(method, avg_slope) %>% dplyr::rename(slope = avg_slope) %>% dplyr::mutate(regression = "LOESS")
) %>% dplyr::mutate(Variable = "Covariates") %>% tidyr::pivot_wider(values_from = slope, names_from = regression, names_prefix = "Slope ") 

p2_covariates <- res %>%
  ggplot(aes(x = n_covariates, y = runtime, color = method)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  scale_y_continuous(transform = "log10") +
  scale_color_manual(values = MY_PALETTE, name = "Method") +
  labs(
    x = "Number of Covariates",
    y = "Runtime (seconds, log scale)"
  ) +
  theme_bw()

# Figure 3: Speedup ####
speedup_data <- res %>%
  group_by(idx) %>%
  mutate(
    base_runtime = runtime[method == BASE_METHOD],
    speedup = base_runtime / runtime
  ) %>%
  ungroup() %>%
  filter(!is.na(speedup), !is.infinite(speedup))

p3_speedup_patients <- speedup_data %>%
  ggplot(aes(x = n_patients, y = speedup, color = method)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_y_log10(
    breaks = c(0.1, 0.5, 1, 2, 5, 10, 50),
    labels = c("0.1×", "0.5×", "1×", "2×", "5×", "10×", "50×")
  ) +
  scale_color_manual(values = MY_PALETTE, name = "Method") +
  scale_linetype_discrete(name = "Number of\nCovariates") +
  labs(
    x = "Number of Patients",
    y = "Speedup relative to NEBULA",
  ) +
  theme_bw()

p3_speedup_covariates <- speedup_data %>%
  ggplot(aes(x = n_covariates, y = speedup, color = method)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  scale_y_log10(
    breaks = c(0.1, 0.5, 1, 2, 5, 10, 50),
    labels = c("0.1×", "0.5×", "1×", "2×", "5×", "10×", "50×")
  ) +
  scale_color_manual(values = MY_PALETTE, name = "Method") +
  scale_linetype_discrete(name = "Number of\nPatients") +
  labs(
    x = "Number of Covariates",
    y = "Speedup relative to NEBULA"
  ) +
  theme_bw()

# Save RDS ####
p1_patients
p2_covariates
p3_speedup_patients
p3_speedup_covariates

saveRDS(p1_patients, "figures/plot_runtime_patients.RDS")
saveRDS(p2_covariates, "figures/plot_runtime_covariates.RDS")
saveRDS(p3_speedup_patients, "figures/plot_speedup_patients.RDS")
saveRDS(p3_speedup_covariates, "figures/plot_speedup_covariates.RDS")

# Supplementary Table: Summary Statistics ####

summary_run_time_patients <- res %>%
  dplyr::group_by(method, n_patients) %>%
  dplyr::summarise(
    mean_runtime = mean(runtime, na.rm = TRUE),
    sd_runtime   = sd(runtime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_runtime = round(mean_runtime, 2),
    sd_runtime   = round(sd_runtime, 2),
    runtime = str_glue("{mean_runtime} ± {sd_runtime}")
  ) %>%
  dplyr::select(method, n_patients, runtime) %>%
  pivot_wider(
    names_from = n_patients,
    values_from = runtime,
    names_prefix = "patients_"
  ) %>%
  dplyr::mutate(method = factor(method, levels = c("devil", "glmGamPoi", "NEBULA"))) %>% 
  dplyr::arrange(method)

summary_run_time_covariates <- res %>%
  dplyr::group_by(method, n_covariates) %>%
  dplyr::summarise(
    mean_runtime = mean(runtime, na.rm = TRUE),
    sd_runtime   = sd(runtime, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mean_runtime = round(mean_runtime, 2),
    sd_runtime   = round(sd_runtime, 2),
    runtime = str_glue("{mean_runtime} ± {sd_runtime}")
  ) %>%
  dplyr::select(method, n_covariates, runtime) %>%
  pivot_wider(
    names_from = n_covariates,
    values_from = runtime,
    names_prefix = "Covariates_"
  ) %>%
  dplyr::mutate(method = factor(method, levels = c("devil", "glmGamPoi", "NEBULA"))) %>% 
  dplyr::arrange(method)

write_csv(summary_run_time_covariates, "summarized_results/supplementary_table_runtime_covariates.csv")
write_csv(summary_run_time_patients, "summarized_results/supplementary_table_runtime_patients.csv")

# Add slope table
write_csv(dplyr::bind_rows(df_slopes_patients, df_slopes_covariates), "summarized_results/supplementary_table_slopes_patients.csv")

