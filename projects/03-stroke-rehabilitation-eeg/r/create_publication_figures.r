#!/usr/bin/env Rscript

# Optional publication-style figures based on the CSV files exported by Python.
# Usage: Rscript r/create_publication_figures.r results

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) args[[1]] else "results"
tables_dir <- file.path(output_dir, "tables")
diagnostics_dir <- file.path(output_dir, "diagnostics")
fig_dir <- file.path(output_dir, "publication_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

packages <- c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "scales", "patchwork", "viridis")
installed <- rownames(installed.packages())
for (pkg in packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
  library(patchwork)
  library(viridis)
})

theme_publication <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size, color = "grey25"),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey95", color = "grey40", linewidth = 0.3),
      strip.text = element_text(face = "bold"),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
      plot.margin = margin(8, 8, 8, 8)
    )
}

save_pub_plot <- function(plot_obj, filename, width = 8, height = 5) {
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot_obj, width = width, height = height)
  ggsave(file.path(fig_dir, paste0(filename, ".png")), plot_obj, width = width, height = height, dpi = 600)
}

plot_confusion_grid <- function(path, title, subtitle, filename) {
  data <- read_csv(path, show_col_types = FALSE) %>%
    mutate(
      stage = factor(stage, levels = c("PRE", "POST")),
      true_label = factor(true_label, levels = c("Left imagery (+1)", "Right imagery (-1)")),
      predicted_label = factor(predicted_label, levels = c("Left imagery (+1)", "Right imagery (-1)")),
      label = as.character(n)
    )

  plot <- ggplot(data, aes(x = predicted_label, y = true_label, fill = n)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = label), size = 4.2, fontface = "bold") +
    facet_grid(stage ~ subject) +
    scale_fill_gradient(low = "#DEEBF7", high = "#08306B", name = "Trials") +
    labs(title = title, subtitle = subtitle, x = "Predicted class", y = "True class") +
    theme_publication(11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
  save_pub_plot(plot, filename, width = 10, height = 6)
}

plot_confusion_grid(
  file.path(tables_dir, "fixed_csp_lda_confusion_matrices_long.csv"),
  "Fixed CSP+LDA Confusion Matrices Across Participants",
  "Held-out TEST trials; rows are true labels and columns are predicted labels",
  "fixed_csp_lda_confusion_matrices"
)

plot_confusion_grid(
  file.path(tables_dir, "training_cv_selected_confusion_matrices_long.csv"),
  "Training-CV Selected Model Confusion Matrices",
  "Model selection used TRAINING data only; final performance uses held-out TEST trials",
  "training_cv_selected_confusion_matrices"
)

trigger_df <- read_csv(file.path(diagnostics_dir, "trigger_onset_diagnostics.csv"), show_col_types = FALSE) %>%
  mutate(
    stage = factor(stage, levels = c("PRE", "POST")),
    run_type = factor(run_type, levels = c("TRAINING", "TEST")),
    imagery_class = factor(imagery_class, levels = c("Left imagery (+1)", "Right imagery (-1)"))
  )

p_trigger_timeline <- ggplot(trigger_df, aes(x = onset_time_s, y = imagery_class, color = imagery_class)) +
  geom_point(size = 1.4, alpha = 0.85) +
  facet_grid(run_type + stage ~ subject, scales = "free_x") +
  labs(
    title = "Trigger Onset Timeline for Left and Right Motor Imagery Trials",
    subtitle = "Each point represents one detected trial onset in the trigger channel",
    x = "Time from file start (s)", y = "Trigger class", color = "Class"
  ) +
  theme_publication(10)
save_pub_plot(p_trigger_timeline, "trigger_onset_timeline", width = 11, height = 7)

trigger_balance <- trigger_df %>% count(subject, stage, run_type, imagery_class, name = "n_trials")
p_trigger_balance <- ggplot(trigger_balance, aes(x = imagery_class, y = n_trials, fill = imagery_class)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.25) +
  geom_text(aes(label = n_trials), vjust = -0.35, size = 3.5) +
  facet_grid(run_type + stage ~ subject) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Class Balance of Trigger Labels Across EEG Files",
    subtitle = "Balanced left/right trial counts reduce the risk of biased classification accuracy",
    x = "Motor imagery class", y = "Number of trials", fill = "Class"
  ) +
  theme_publication(10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
save_pub_plot(p_trigger_balance, "trigger_class_balance", width = 11, height = 7)

electrode_activity <- read_csv(
  file.path(diagnostics_dir, "electrode_activity_summary_fixed_window.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    stage = factor(stage, levels = c("PRE", "POST")),
    run_type = factor(run_type, levels = c("TRAINING", "TEST")),
    imagery_class = factor(imagery_class, levels = c("Left imagery (+1)", "Right imagery (-1)")),
    channel = factor(channel, levels = unique(channel[order(channel_index)]))
  )

electrode_test <- electrode_activity %>% filter(run_type == "TEST")
p_electrode_heatmap <- ggplot(electrode_test, aes(x = channel, y = subject, fill = mean_log_variance)) +
  geom_tile(color = "white", linewidth = 0.35) +
  facet_grid(imagery_class ~ stage) +
  scale_fill_viridis_c(option = "C", name = "Mean log\nvariance") +
  labs(
    title = "Electrode-Level Activity During Motor Imagery TEST Trials",
    subtitle = "Mean log-variance in the fixed 2–8 s, 8–30 Hz analysis window",
    x = "EEG electrode", y = "Participant"
  ) +
  theme_publication(10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
save_pub_plot(p_electrode_heatmap, "electrode_activity_heatmap", width = 11, height = 6)

electrode_rank <- electrode_test %>%
  group_by(subject, stage, channel) %>%
  summarise(mean_log_variance = mean(mean_log_variance), .groups = "drop") %>%
  group_by(subject, stage) %>%
  slice_max(mean_log_variance, n = 6, with_ties = FALSE) %>%
  ungroup()

p_top_electrodes <- ggplot(electrode_rank, aes(x = reorder(channel, mean_log_variance), y = mean_log_variance)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  coord_flip() +
  facet_grid(stage ~ subject, scales = "free_y") +
  labs(
    title = "Highest-Activity Electrodes During TEST Motor Imagery",
    subtitle = "Top six electrodes per participant and session using mean log-variance",
    x = "EEG electrode", y = "Mean log-variance"
  ) +
  theme_publication(10) +
  theme(legend.position = "none")
save_pub_plot(p_top_electrodes, "top_electrodes_by_mean_log_variance", width = 10, height = 7)

fixed_results <- read_csv(file.path(tables_dir, "stroke_bci_fixed_csp_lda_results.csv"), show_col_types = FALSE) %>%
  mutate(stage = factor(stage, levels = c("pre", "post"), labels = c("PRE", "POST")))

p_fixed_accuracy <- ggplot(fixed_results, aes(x = stage, y = test_accuracy_percent, group = subject)) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.45, color = "grey35") +
  geom_line(aes(color = subject), linewidth = 0.9) +
  geom_point(aes(color = subject), size = 3.2) +
  scale_y_continuous(limits = c(45, 100), breaks = seq(50, 100, 10)) +
  labs(
    title = "Fixed CSP+LDA Motor Imagery Decoding Accuracy",
    subtitle = "Same 2–8 s window, 8–30 Hz band, CSP(4), and LDA classifier for all sessions",
    x = "Rehabilitation session", y = "Held-out TEST accuracy (%)", color = "Participant"
  ) +
  theme_publication(12)
save_pub_plot(p_fixed_accuracy, "fixed_csp_lda_accuracy_pre_post", width = 7.5, height = 5)

message("Publication figures saved to: ", normalizePath(fig_dir, mustWork = FALSE))
