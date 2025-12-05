
rm(list=ls())
library(tidyverse)
library(ggsci)

format_human <- function(x) {
  ifelse(abs(x) >= 1e9,  paste0(round(x / 1e9, 1), "B"),
         ifelse(abs(x) >= 1e6,  paste0(round(x / 1e6, 1), "M"),
                ifelse(abs(x) >= 1e3,  paste0(round(x / 1e3, 1), "k"),
                       as.character(x))))
}

df = read.delim("results/memory_gpu.csv", sep = ",") %>%
  dplyr::mutate(memory_gb = 0.0009765625 * gpu_used_memory_mb)

df$n_genes = factor(format_human(df$n_genes), levels = c("1k", "5k", "10k"))
df$n_cells = factor(format_human(df$n_cells), levels = c("10k", "100k", "1M"))

# --- summarise memory usage by groups (needed for error bars) ---
df_summary <- df %>%
  group_by(n_genes, n_cells, batch_size) %>%
  summarise(
    mean_mem = mean(memory_gb, na.rm = TRUE),
    sd_mem   = ifelse(dplyr::n() > 1, sd(memory_gb, na.rm = TRUE), 0),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_mem - sd_mem,
    ymax = mean_mem + sd_mem
  )

# --- plot ---
plot_memory_usage = ggplot(df_summary,
       aes(x = factor(batch_size),
           y = mean_mem, fill = mean_mem)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax),
                position = position_dodge(width = 0.9),
                width = 0.2) +
  ggh4x::facet_nested("N genes" + n_genes ~ "N cells" + n_cells, scales = "free_y") +
  labs(
    x = "Batch size",
    y = "Used GPU memory (GB)"
  ) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_fill_gradient(low = "forestgreen", high = "firebrick")

ggsave(plot = plot_memory_usage, filename = "figures/memory_usage.pdf", width = 6, height = 6)
saveRDS(plot_memory_usage, "figures/memory_usage.RDS")

# Strong Scaling ####
df = read.delim("results/scaling_gpu.csv", sep = ",")
df$n_genes = factor(format_human(df$n_genes), levels = c("1k", "5k", "10k"))
df$n_cells = factor(format_human(df$n_cells), levels = c("10k", "100k", "1M"))

plot_strong_scaling = ggplot(df,
       aes(x = n_gpus,
           y = speedup, col = as.factor(batch_size))) +
  geom_point() +
  geom_line() +
  ggh4x::facet_nested("N genes" + n_genes ~ "N cells" + n_cells) +
  labs(
    x = "Number of GPUs",
    y = "Speedup", colour = "Batch size"
  ) +
  theme_bw() +
  ggsci::scale_color_bmj()


ggsave(plot = plot_strong_scaling, filename = "figures/strong_scaling.pdf", width = 6, height = 6)
saveRDS(plot_strong_scaling, "figures/strong_scaling.RDS")

# Both

pboth = plot_strong_scaling + plot_memory_usage + patchwork::plot_layout(design = "A\nB") +
  patchwork::plot_annotation(tag_levels = c("A"))


ggsave(plot = pboth, filename = "figures/gpu_memory_and_strongscaling.pdf", width = 6, height = 10)
saveRDS(pboth, "figures/gpu_memory_and_strongscaling.RDS")
