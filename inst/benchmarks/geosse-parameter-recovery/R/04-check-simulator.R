# Check exact-tip GeoSSE simulations without fitting the models.
#
# Usage:
# Rscript inst/benchmarks/geosse-parameter-recovery/R/04-check-simulator.R 3

library(divextra)
library(diversitree)
library(ape)

source(file.path(
  "inst", "benchmarks", "geosse-parameter-recovery",
  "R", "benchmark-functions.R"
))

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) == 0) {
  pilot.replicates <- 3L
} else {
  pilot.replicates <- as.integer(arguments[1])
}

if (is.na(pilot.replicates) || pilot.replicates < 1) {
  stop("The number of pilot replicates must be a positive integer")
}

scenarios <- read_benchmark_scenarios(calibrated = TRUE)
pilot.rows <- list()
row.number <- 0L

for (scenario.number in seq_len(nrow(scenarios))) {
  scenario <- scenarios[scenario.number, , drop = FALSE]

  for (replicate.number in seq_len(pilot.replicates)) {
    message(
      "Checking ", scenario$scenario_id,
      ", replicate ", replicate.number
    )

    start.time <- proc.time()["elapsed"]

    simulation <- simulate_benchmark_dataset(
      scenario = scenario,
      replicate.number = replicate.number
    )

    elapsed.seconds <- unname(proc.time()["elapsed"] - start.time)
    terminal.lineages <- !simulation$table$split
    extinct.terminals <- terminal.lineages & simulation$table$extinct

    row.number <- row.number + 1L
    pilot.rows[[row.number]] <- data.frame(
      scenario_id = scenario$scenario_id,
      replicate = replicate.number,
      target_tips = scenario$target_tips,
      observed_tips = ape::Ntip(simulation$tree),
      exact_tip_count = ape::Ntip(simulation$tree) == scenario$target_tips,
      attempts = simulation$attempt,
      total_extinction_rejections =
        unname(simulation$rejection.counts["total_extinction"]),
      epoch1_target_rejections =
        unname(simulation$rejection.counts["target_reached_in_epoch1"]),
      other_rejections = sum(simulation$rejection.counts) -
        unname(simulation$rejection.counts["total_extinction"]) -
        unname(simulation$rejection.counts["target_reached_in_epoch1"]),
      extinct_terminal_lineages = sum(extinct.terminals),
      epoch1_duration = simulation$epoch1.duration,
      simulated_time = simulation$simulated.time,
      tree_height = simulation$tree.height,
      boundary_age = simulation$boundary.age,
      realized_boundary_fraction = simulation$realized.boundary.fraction,
      elapsed_seconds = elapsed.seconds,
      stringsAsFactors = FALSE
    )
  }
}

pilot <- do.call(rbind, pilot.rows)
rownames(pilot) <- NULL

output.dir <- file.path(
  benchmark_directory(), "results", "simulator-check"
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  pilot,
  file.path(output.dir, "simulator-pilot.csv"),
  row.names = FALSE
)

scenario.ids <- unique(pilot$scenario_id)
summary.rows <- vector("list", length(scenario.ids))

for (scenario.number in seq_along(scenario.ids)) {
  scenario.id <- scenario.ids[scenario.number]
  scenario.data <- pilot[pilot$scenario_id == scenario.id, , drop = FALSE]

  summary.rows[[scenario.number]] <- data.frame(
    scenario_id = scenario.id,
    n_simulations = nrow(scenario.data),
    target_tips = scenario.data$target_tips[1],
    all_tip_counts_exact = all(scenario.data$exact_tip_count),
    median_attempts = median(scenario.data$attempts),
    median_extinct_lineages = median(
      scenario.data$extinct_terminal_lineages
    ),
    median_realized_boundary_fraction = median(
      scenario.data$realized_boundary_fraction
    ),
    median_runtime_seconds = median(scenario.data$elapsed_seconds),
    stringsAsFactors = FALSE
  )
}

pilot.summary <- do.call(rbind, summary.rows)
rownames(pilot.summary) <- NULL

write.csv(
  pilot.summary,
  file.path(output.dir, "simulator-pilot-summary.csv"),
  row.names = FALSE
)

report.file <- file.path(output.dir, "simulator-pilot-report.md")
table.lines <- c(
  "| Scenario | N simulations | Target tips | Exact | Median attempts | Median extinct | Median boundary fraction | Median runtime (s) |",
  "|---|---:|---:|:---:|---:|---:|---:|---:|"
)

for (scenario.number in seq_len(nrow(pilot.summary))) {
  result <- pilot.summary[scenario.number, , drop = FALSE]

  table.lines <- c(
    table.lines,
    paste0(
      "| ", result$scenario_id,
      " | ", result$n_simulations,
      " | ", result$target_tips,
      " | ", result$all_tip_counts_exact,
      " | ", format(result$median_attempts, trim = TRUE),
      " | ", format(result$median_extinct_lineages, trim = TRUE),
      " | ", sprintf("%.3f", result$median_realized_boundary_fraction),
      " | ", sprintf("%.3f", result$median_runtime_seconds),
      " |"
    )
  )
}

report.lines <- c(
  "# Exact-tip GeoSSE simulator pilot",
  "",
  paste0("Pilot replicates per scenario: ", pilot.replicates),
  "",
  paste0("All accepted trees had exact target size: ",
         all(pilot$exact_tip_count)),
  "",
  "The CSV files in this directory contain the complete numerical results.",
  "The realized boundary fraction is measured on the reconstructed crown tree.",
  "",
  "## Scenario summary",
  "",
  table.lines,
  "",
  "## Interpretation",
  "",
  paste0(
    "- Exact target size was recovered in ", sum(pilot$exact_tip_count),
    " of ", nrow(pilot), " accepted simulations."
  ),
  paste0(
    "- Accepted trees contained ", min(pilot$extinct_terminal_lineages),
    " to ", max(pilot$extinct_terminal_lineages),
    " extinct terminal lineages."
  ),
  paste0(
    "- Accepted simulations required ", min(pilot$attempts),
    " to ", max(pilot$attempts), " attempts."
  ),
  "- The realized epoch fraction varies because final tree age is stochastic.",
  "  Boundary calibration should be assessed separately from exact tip count."
)

writeLines(report.lines, report.file)

print(pilot.summary, row.names = FALSE)
message("Saved simulator checks in: ", output.dir)
