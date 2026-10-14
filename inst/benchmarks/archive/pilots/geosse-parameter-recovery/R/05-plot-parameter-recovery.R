# Plot episodic GeoSSE parameter recovery.
# Run after 03-summarize.R.

library(ggplot2)

benchmark.dir <- file.path(
  "inst", "benchmarks", "archive", "pilots", "geosse-parameter-recovery"
)
summary.dir <- file.path(benchmark.dir, "results", "summaries")
plot.dir <- file.path(benchmark.dir, "results", "plots")
dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)

datasets <- read.csv(
  file.path(summary.dir, "dataset-estimates.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
recovery <- read.csv(
  file.path(summary.dir, "parameter-recovery.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

parameter.names <- c(
  "s_recent", "d_recent", "x_recent",
  "s_old", "d_old", "x_old"
)
parameter.labels <- c(
  s_recent = "Speciation, recent",
  d_recent = "Dispersal, recent",
  x_recent = "Contraction, recent",
  s_old = "Speciation, old",
  d_old = "Dispersal, old",
  x_old = "Contraction, old"
)
scenario.order <- c(
  "P01", "P02", "P03", "P04", "P05",
  "E25", "E75", "N050", "N500"
)

point.rows <- list()
row.number <- 0L

for (parameter in parameter.names) {
  true.column <- paste0("true_", parameter)
  estimate.column <- paste0("estimate_", parameter)

  for (dataset.number in seq_len(nrow(datasets))) {
    row.number <- row.number + 1L
    point.rows[[row.number]] <- data.frame(
      scenario_id = datasets$scenario_id[dataset.number],
      replicate = datasets$replicate[dataset.number],
      parameter = parameter,
      truth = datasets[[true.column]][dataset.number],
      estimate = datasets[[estimate.column]][dataset.number],
      stringsAsFactors = FALSE
    )
  }
}

points <- do.call(rbind, point.rows)
rownames(points) <- NULL

points$scenario_id <- factor(points$scenario_id, levels = scenario.order)
points$parameter <- factor(
  points$parameter,
  levels = parameter.names,
  labels = parameter.labels[parameter.names]
)
points$recovery_ratio <- points$estimate / points$truth
recovery$scenario_id <- factor(recovery$scenario_id, levels = scenario.order)
recovery$parameter <- factor(
  recovery$parameter,
  levels = parameter.names,
  labels = parameter.labels[parameter.names]
)

scenario.colors <- c(
  P01 = "#000000", P02 = "#0072B2", P03 = "#009E73",
  P04 = "#D55E00", P05 = "#CC79A7", E25 = "#56B4E9",
  E75 = "#E69F00", N050 = "#999999", N500 = "#332288"
)

figure.theme <- theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text = element_text(color = "black")
  )

recovery.scatter <- ggplot(
  points,
  aes(x = scenario_id, y = recovery_ratio, color = scenario_id)
) +
  geom_hline(
    yintercept = 1, linewidth = 0.4,
    linetype = "dashed", color = "grey40"
  ) +
  geom_point(
    size = 2,
    alpha = 0.85,
    position = position_jitter(width = 0.12, height = 0, seed = 1)
  ) +
  facet_wrap(~parameter, ncol = 3) +
  scale_color_manual(values = scenario.colors, drop = FALSE) +
  scale_y_log10() +
  labs(
    x = "Simulation scenario",
    y = "Estimate / simulated value",
    caption = paste0(
      "Local pilot only: ", nrow(datasets),
      " fitted datasets; logarithmic axis; dashed line indicates perfect recovery."
    )
  ) +
  figure.theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

ggsave(
  file.path(plot.dir, "parameter-recovery-estimated-vs-true.pdf"),
  recovery.scatter, width = 7.2, height = 5.2, units = "in"
)
ggsave(
  file.path(plot.dir, "parameter-recovery-estimated-vs-true.png"),
  recovery.scatter, width = 7.2, height = 5.2, units = "in", dpi = 300
)

recovery.intervals <- ggplot(recovery, aes(y = scenario_id)) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    width = 0, linewidth = 0.55, color = "grey45",
    orientation = "y"
  ) +
  geom_point(aes(x = mean, color = scenario_id), size = 2.2) +
  geom_point(
    aes(x = truth), shape = 3, size = 2.8,
    stroke = 0.8, color = "black"
  ) +
  facet_wrap(~parameter, ncol = 3, scales = "free_x") +
  scale_x_log10() +
  scale_color_manual(values = scenario.colors, guide = "none") +
  labs(
    x = "Estimated rate",
    y = NULL,
    caption = paste0(
      "Colored points are means; horizontal lines are empirical 95% ranges; ",
      "black crosses are simulated values. Pilot n = 1-3 per scenario."
    )
  ) +
  figure.theme

ggsave(
  file.path(plot.dir, "parameter-recovery-by-scenario.pdf"),
  recovery.intervals, width = 7.2, height = 6.0, units = "in"
)
ggsave(
  file.path(plot.dir, "parameter-recovery-by-scenario.png"),
  recovery.intervals, width = 7.2, height = 6.0, units = "in", dpi = 300
)

message("Saved parameter-recovery plots in: ", plot.dir)
