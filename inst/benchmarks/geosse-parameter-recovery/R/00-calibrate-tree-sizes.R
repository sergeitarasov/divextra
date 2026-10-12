# Calibrate fixed simulation durations for approximate tree-size classes.

library(divextra)
library(diversitree)
library(ape)

source(file.path(
  "inst", "benchmarks", "geosse-parameter-recovery",
  "R", "benchmark-functions.R"
))

benchmark.dir <- benchmark_directory()
scenarios <- read_benchmark_scenarios(calibrated = FALSE)
scenarios$scenario_number <- seq_len(nrow(scenarios))

calibration.dir <- file.path(
  benchmark.dir,
  "results", "calibrated"
)

dir.create(calibration.dir, recursive = TRUE, showWarnings = FALSE)

calibrated.rows <- vector("list", nrow(scenarios))

for (scenario.number in seq_len(nrow(scenarios))) {
  scenario <- scenarios[scenario.number, , drop = FALSE]

  message(
    "Calibrating ", scenario$scenario_id,
    ": target median = ", scenario$target_tips, " tips"
  )

  calibration <- calibrate_scenario_duration(
    scenario = scenario,
    calibration.replicates = 10,
    calibration.steps = 7
  )

  calibrated.rows[[scenario.number]] <- cbind(
    scenario,
    calibration
  )

  message(
    "  total time = ", round(calibration$total_time, 3),
    "; validation median = ", calibration$calibration_median_tips
  )
}

calibrated.scenarios <- do.call(rbind, calibrated.rows)
rownames(calibrated.scenarios) <- NULL

output.file <- file.path(
  calibration.dir,
  "scenarios-calibrated.csv"
)

write.csv(
  calibrated.scenarios,
  output.file,
  row.names = FALSE
)

message("Saved: ", output.file)

