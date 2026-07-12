#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)

df <- fread("/path/to/Desktop/benchmarking/tables/pass_depth.csv")
df_new <- df %>%
  group_by(chromosome, position, target) %>%
  summarise(mean_depth = mean(depth, na.rm = TRUE), .groups = "drop")
nrow(df)
nrow(df_new)

color_map <- c("on_target" = "tomato", "off_target" = "gray")

ggplot(df_new, aes(x = position, y = mean_depth, color = target)) +
  geom_point(alpha = 0.4, size = 1) +  # Adjusted point size for clarity
  scale_color_manual(
    values = color_map,
    name = "Target status",
    labels = c("Off target", "On target")
  ) +
  scale_x_continuous(labels = NULL, breaks = NULL) +  # Remove x-axis labels
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "transparent"),
    axis.title.y = element_text(size = 14, angle = 90),
    axis.title.x = element_text(size = 14),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
    panel.grid.minor.y = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +  # Larger legend points
  labs(
    x = "Genome position",
    y = "Depth"
  )


ggsave("depth_plot_barcode.png", width = 15, height = 5, dpi = 400)