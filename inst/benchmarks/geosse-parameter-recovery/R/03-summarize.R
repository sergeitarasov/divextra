# Summarize completed episodic GeoSSE recovery results.

source(file.path(
  "inst", "benchmarks", "geosse-parameter-recovery",
  "R", "benchmark-functions.R"
))

benchmark.dir <- benchmark_directory()
raw.dir <- file.path(benchmark.dir, "results", "raw")
summary.dir <- file.path(benchmark.dir, "results", "summaries")

dir.create(summary.dir, recursive = TRUE, showWarnings = FALSE)

result.files <- list.files(
  raw.dir,
  pattern = "[.]rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(result.files) == 0) {
  stop("No benchmark result files were found in ", raw.dir)
}

dataset.rows <- vector("list", length(result.files))
start.rows <- list()
start.row.number <- 0L

parameter.names <- c(
  "s_recent", "d_recent", "x_recent",
  "s_old", "d_old", "x_old"
)

for (file.number in seq_along(result.files)) {
  result <- readRDS(result.files[file.number])

  truth <- c(
    s_recent = result$truth["recent", "sA"],
    d_recent = result$truth["recent", "dA"],
    x_recent = result$truth["recent", "xA"],
    s_old = result$truth["old", "sA"],
    d_old = result$truth["old", "dA"],
    x_old = result$truth["old", "xA"]
  )

  estimates <- result$episodic_estimates

  episodic.logLik <- NA_real_
  episodic.AIC <- NA_real_
  constant.logLik <- NA_real_
  constant.AIC <- NA_real_

  if (!is.null(result$episodic$best)) {
    episodic.logLik <- result$episodic$best$lnLik
    episodic.AIC <- result$episodic$best$AIC
  }

  if (!is.null(result$constant$best)) {
    constant.logLik <- result$constant$best$lnLik
    constant.AIC <- result$constant$best$AIC
  }

  episodic.successful <- result$episodic$starts[
    vapply(
      result$episodic$starts,
      function(x) isTRUE(x$success),
      logical(1)
    )
  ]

  if (length(episodic.successful) > 1) {
    start.logLik <- vapply(
      episodic.successful,
      function(x) x$lnLik,
      numeric(1)
    )
    start.range <- max(start.logLik) - min(start.logLik)
  } else {
    start.range <- NA_real_
  }

  dataset.row <- data.frame(
    scenario_id = result$scenario_id,
    replicate = result$replicate,
    seed = result$seed,
    n_tip = result$n_tip,
    target_tips = result$scenario$target_tips,
    epoch1_duration = result$epoch1_duration,
    simulated_time = result$simulated_time,
    tree_height = result$tree_height,
    boundary_age = result$boundary_age,
    nominal_boundary_fraction = result$nominal_boundary_fraction,
    realized_boundary_fraction = result$realized_boundary_fraction,
    simulation_attempt = result$simulation_attempt,
    simulated_root_state = result$simulated_root_state,
    reconstructed_root_state = result$reconstructed_root_state,
    elapsed_seconds = result$elapsed_seconds,
    episodic_logLik = episodic.logLik,
    constant_logLik = constant.logLik,
    episodic_AIC = episodic.AIC,
    constant_AIC = constant.AIC,
    delta_AIC_constant_minus_episodic = constant.AIC - episodic.AIC,
    episodic_best_start = result$episodic$best.start,
    constant_best_start = result$constant$best.start,
    episodic_start_logLik_range = start.range,
    stringsAsFactors = FALSE
  )

  for (parameter in parameter.names) {
    dataset.row[[paste0("true_", parameter)]] <- truth[parameter]
    dataset.row[[paste0("estimate_", parameter)]] <- estimates[parameter]
  }

  dataset.rows[[file.number]] <- dataset.row

  models <- c("episodic", "constant")

  for (model in models) {
    model.result <- result[[model]]

    for (start.number in seq_along(model.result$starts)) {
      fit <- model.result$starts[[start.number]]
      start.row.number <- start.row.number + 1L

      start.rows[[start.row.number]] <- data.frame(
        scenario_id = result$scenario_id,
        replicate = result$replicate,
        model = model,
        start = fit$start,
        eps = fit$eps,
        success = fit$success,
        lnLik = if (isTRUE(fit$success)) fit$lnLik else NA_real_,
        AIC = if (isTRUE(fit$success)) fit$AIC else NA_real_,
        convergence = if (isTRUE(fit$success)) fit$convergence else NA_integer_,
        elapsed_seconds = fit$elapsed_seconds,
        error = fit$error,
        stringsAsFactors = FALSE
      )
    }
  }
}

datasets <- do.call(rbind, dataset.rows)
starts <- do.call(rbind, start.rows)
rownames(datasets) <- NULL
rownames(starts) <- NULL

write.csv(
  datasets,
  file.path(summary.dir, "dataset-estimates.csv"),
  row.names = FALSE
)

write.csv(
  starts,
  file.path(summary.dir, "optimization-starts.csv"),
  row.names = FALSE
)

scenario.ids <- unique(datasets$scenario_id)
recovery.rows <- list()
recovery.row.number <- 0L

for (scenario.id in scenario.ids) {
  scenario.data <- datasets[
    datasets$scenario_id == scenario.id,
    ,
    drop = FALSE
  ]

  for (parameter in parameter.names) {
    true.column <- paste0("true_", parameter)
    estimate.column <- paste0("estimate_", parameter)

    true.value <- scenario.data[[true.column]][1]
    estimates <- scenario.data[[estimate.column]]
    estimates <- estimates[is.finite(estimates)]

    recovery.row.number <- recovery.row.number + 1L

    if (length(estimates) == 0) {
      recovery.rows[[recovery.row.number]] <- data.frame(
        scenario_id = scenario.id,
        parameter = parameter,
        truth = true.value,
        n = 0,
        mean = NA_real_,
        median = NA_real_,
        bias = NA_real_,
        relative_bias = NA_real_,
        RMSE = NA_real_,
        lower_95 = NA_real_,
        upper_95 = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    recovery.rows[[recovery.row.number]] <- data.frame(
      scenario_id = scenario.id,
      parameter = parameter,
      truth = true.value,
      n = length(estimates),
      mean = mean(estimates),
      median = median(estimates),
      bias = mean(estimates - true.value),
      relative_bias = mean((estimates - true.value) / true.value),
      RMSE = sqrt(mean((estimates - true.value)^2)),
      lower_95 = unname(quantile(estimates, 0.025, type = 7)),
      upper_95 = unname(quantile(estimates, 0.975, type = 7)),
      stringsAsFactors = FALSE
    )
  }
}

recovery <- do.call(rbind, recovery.rows)
rownames(recovery) <- NULL

write.csv(
  recovery,
  file.path(summary.dir, "parameter-recovery.csv"),
  row.names = FALSE
)

scenario.summary <- data.frame(
  scenario_id = scenario.ids,
  n_datasets = NA_integer_,
  median_n_tip = NA_real_,
  median_runtime_minutes = NA_real_,
  episodic_fit_success = NA_real_,
  constant_fit_success = NA_real_,
  starts_agree = NA_real_,
  episodic_AIC_support = NA_real_,
  stringsAsFactors = FALSE
)

for (scenario.number in seq_along(scenario.ids)) {
  scenario.id <- scenario.ids[scenario.number]
  scenario.data <- datasets[
    datasets$scenario_id == scenario.id,
    ,
    drop = FALSE
  ]

  scenario.summary$n_datasets[scenario.number] <- nrow(scenario.data)
  scenario.summary$median_n_tip[scenario.number] <- median(
    scenario.data$n_tip
  )
  scenario.summary$median_runtime_minutes[scenario.number] <- median(
    scenario.data$elapsed_seconds / 60
  )
  scenario.summary$episodic_fit_success[scenario.number] <- mean(
    is.finite(scenario.data$episodic_logLik)
  )
  scenario.summary$constant_fit_success[scenario.number] <- mean(
    is.finite(scenario.data$constant_logLik)
  )
  scenario.summary$starts_agree[scenario.number] <- mean(
    scenario.data$episodic_start_logLik_range < 1e-4,
    na.rm = TRUE
  )
  scenario.summary$episodic_AIC_support[scenario.number] <- mean(
    scenario.data$delta_AIC_constant_minus_episodic > 2,
    na.rm = TRUE
  )
}

write.csv(
  scenario.summary,
  file.path(summary.dir, "scenario-summary.csv"),
  row.names = FALSE
)

message("Saved summaries in: ", summary.dir)

