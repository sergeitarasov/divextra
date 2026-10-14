# Run a resumable local batch.
#
# Usage:
# Rscript inst/benchmarks/archive/pilots/geosse-parameter-recovery/R/02-run-local.R 5

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) > 1) {
  stop("Supply only the number of replicates per scenario")
}

if (length(arguments) == 0) {
  replicates.per.scenario <- 1L
} else {
  replicates.per.scenario <- as.integer(arguments[1])
}

if (is.na(replicates.per.scenario) || replicates.per.scenario < 1) {
  stop("The number of replicates must be a positive integer")
}

benchmark.dir <- file.path(
  "inst", "benchmarks", "archive", "pilots", "geosse-parameter-recovery"
)

scenario.file <- file.path(
  benchmark.dir,
  "results", "calibrated", "scenarios-calibrated.csv"
)

if (!file.exists(scenario.file)) {
  stop(
    "Run 00-calibrate-tree-sizes.R before running the local benchmark"
  )
}

scenarios <- read.csv(
  scenario.file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

worker.script <- file.path(
  benchmark.dir,
  "R", "01-run-one.R"
)

for (scenario.number in seq_len(nrow(scenarios))) {
  scenario.id <- scenarios$scenario_id[scenario.number]

  for (replicate.number in seq_len(replicates.per.scenario)) {
    message(
      "Starting ", scenario.id,
      ", replicate ", replicate.number
    )

    status <- system2(
      "Rscript",
      args = c(worker.script, scenario.id, replicate.number)
    )

    if (status != 0) {
      warning(
        "Failed: ", scenario.id,
        ", replicate ", replicate.number
      )
    }
  }
}

