suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(here)
  library(scales)
})

table_dir <- here::here("results", "tables")
output_dir <- Sys.getenv(
  "NOSOCOMIAL_FIGURE_DIR",
  here::here("results", "figures")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

baseline <- read.csv(file.path(table_dir, "testing-rate-r0-1-03.csv")) |>
  mutate(scenario = "R0 = 1.03")

high_r0 <- read.csv(file.path(table_dir, "testing-rate-r0-2-5.csv")) |>
  mutate(scenario = "R0 = 2.5")

results <- bind_rows(baseline, high_r0) |>
  mutate(scenario = factor(scenario, levels = c("R0 = 1.03", "R0 = 2.5")))

means <- results |>
  group_by(scenario, test_rate) |>
  summarise(
    se = if (n() > 1L) sd(mode_acc) / sqrt(n()) else 0,
    mode_acc = mean(mode_acc),
    .groups = "drop"
  )

palette <- c("R0 = 1.03" = "#3498db", "R0 = 2.5" = "#e67e22")

p <- ggplot() +
  geom_jitter(
    data = results,
    aes(test_rate, mode_acc, colour = scenario),
    width = 0.004,
    height = 0,
    alpha = 0.10,
    size = 0.8
  ) +
  geom_ribbon(
    data = means,
    aes(test_rate, ymin = mode_acc - se, ymax = mode_acc + se,
        fill = scenario),
    alpha = 0.16,
    colour = NA
  ) +
  geom_smooth(
    data = results,
    aes(test_rate, mode_acc, colour = scenario),
    method = "loess",
    span = 0.35,
    se = FALSE,
    linewidth = 1.2
  ) +
  geom_point(
    data = means,
    aes(test_rate, mode_acc, colour = scenario),
    size = 2.4
  ) +
  scale_colour_manual(
    values = palette,
    labels = c(expression(R[0] == 1.03), expression(R[0] == 2.5)),
    name = NULL
  ) +
  scale_fill_manual(values = palette, guide = "none") +
  scale_x_continuous(
    breaks = seq(0, 1, by = 0.2),
    labels = percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = percent_format(accuracy = 1)
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = expression(p[d]), y = "MAP accuracy") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.88, 0.12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#e5e5e5"),
    axis.line = element_line(colour = "#cccccc"),
    axis.ticks = element_line(colour = "#cccccc")
  )

ggsave(
  file.path(output_dir, "testing-rate-vs-accuracy.png"),
  p,
  width = 10,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(output_dir, "testing-rate-vs-accuracy.pdf"),
  p,
  width = 10,
  height = 6
)
