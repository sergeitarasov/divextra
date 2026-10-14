# Simulate and fit one episodic GeoSSE benchmark dataset.
#
# Usage:
# Rscript inst/benchmarks/archive/pilots/geosse-parameter-recovery/R/01-run-one.R P02 1

library(divextra)
library(diversitree)
library(ape)

source(file.path(
  "inst", "benchmarks", "archive", "pilots", "geosse-parameter-recovery",
  "R", "benchmark-functions.R"
))

arguments <- commandArgs(trailingOnly = TRUE)

if (length(arguments) != 2) {
  stop("Supply scenario_id and replicate number")
}

scenario.id <- arguments[1]
replicate.number <- as.integer(arguments[2])

if (is.na(replicate.number) || replicate.number < 1) {
  stop("Replicate number must be a positive integer")
}

benchmark.dir <- benchmark_directory()
scenarios <- read_benchmark_scenarios(calibrated = TRUE)

scenario.row <- which(scenarios$scenario_id == scenario.id)

if (length(scenario.row) != 1) {
  stop("Unknown or duplicated scenario_id: ", scenario.id)
}

scenario <- scenarios[scenario.row, , drop = FALSE]

output.dir <- file.path(
  benchmark.dir,
  "results", "raw", scenario.id
)

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

output.file <- file.path(
  output.dir,
  sprintf("replicate-%03d.rds", replicate.number)
)

if (file.exists(output.file)) {
  message("Result already exists; skipping: ", output.file)
  quit(save = "no", status = 0)
}

message(
  "Scenario ", scenario.id,
  ", replicate ", replicate.number
)

analysis.start <- proc.time()["elapsed"]

simulation <- simulate_benchmark_dataset(
  scenario = scenario,
  replicate.number = replicate.number
)

tree <- simulation$tree

message(
  "  simulated ", ape::Ntip(tree),
  " extant tips after ", simulation$attempt,
  " attempt(s)"
)

episodic.model.file <- file.path(
  benchmark.dir,
  "config", "model-episodic.yml"
)

constant.model.file <- file.path(
  benchmark.dir,
  "config", "model-constant.yml"
)

episodic.fit <- fit_benchmark_model(
  tree = tree,
  boundary.age = simulation$boundary.age,
  model.file = episodic.model.file,
  model.name = "episodic",
  eps.values = c(0.1, 0.4, 0.8)
)

constant.fit <- fit_benchmark_model(
  tree = tree,
  boundary.age = simulation$boundary.age,
  model.file = constant.model.file,
  model.name = "constant",
  eps.values = c(0.1, 0.4, 0.8)
)

git.commit <- tryCatch(
  system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
  error = function(error) NA_character_
)

if (length(git.commit) != 1) {
  git.commit <- NA_character_
}

analysis.elapsed <- unname(proc.time()["elapsed"] - analysis.start)

result <- list(
  scenario = scenario,
  scenario_id = scenario.id,
  replicate = replicate.number,
  seed = simulation$seed,
  simulation_attempt = simulation$attempt,
  truth = simulation$parameters$geosse,
  simulation_table = simulation$table,
  tree = tree,
  n_tip = ape::Ntip(tree),
  epoch1_duration = simulation$epoch1.duration,
  simulated_time = simulation$simulated.time,
  tree_height = simulation$tree.height,
  boundary_age = simulation$boundary.age,
  nominal_boundary_fraction = scenario$boundary_fraction,
  realized_boundary_fraction = simulation$realized.boundary.fraction,
  simulated_root_state = simulation$simulated.root.state,
  reconstructed_root_state = simulation$reconstructed.root.state,
  episodic = episodic.fit,
  constant = constant.fit,
  episodic_estimates = extract_epoch_estimates(episodic.fit),
  constant_estimates = extract_epoch_estimates(constant.fit),
  elapsed_seconds = analysis.elapsed,
  git_commit = git.commit,
  versions = package_version_record()
)

temporary.file <- paste0(output.file, ".tmp")
saveRDS(result, temporary.file, compress = "xz")
file.rename(temporary.file, output.file)

message(
  "  saved in ", round(analysis.elapsed / 60, 2),
  " minutes: ", output.file
)

