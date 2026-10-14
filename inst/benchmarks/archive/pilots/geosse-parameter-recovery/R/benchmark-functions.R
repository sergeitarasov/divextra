# Shared functions for the episodic GeoSSE parameter-recovery benchmark.

benchmark_directory <- function() {
  file.path("inst", "benchmarks", "archive", "pilots", "geosse-parameter-recovery")
}

read_benchmark_scenarios <- function(calibrated = TRUE) {
  benchmark.dir <- benchmark_directory()

  if (calibrated) {
    scenario.file <- file.path(
      benchmark.dir,
      "results", "calibrated", "scenarios-calibrated.csv"
    )
  } else {
    scenario.file <- file.path(benchmark.dir, "config", "scenarios.csv")
  }

  if (!file.exists(scenario.file)) {
    stop("Scenario file does not exist: ", scenario.file)
  }

  read.csv(scenario.file, stringsAsFactors = FALSE, check.names = FALSE)
}

make_simulation_parameters <- function(scenario) {
  parameter.names <- diversitree:::default.argnames.geosse()

  old.parameters <- c(
    sA = scenario$s_old,
    sB = scenario$s_old,
    sAB = scenario$s_old,
    xA = scenario$x_old,
    xB = scenario$x_old,
    dA = scenario$d_old,
    dB = scenario$d_old
  )

  recent.parameters <- c(
    sA = scenario$s_recent,
    sB = scenario$s_recent,
    sAB = scenario$s_recent,
    xA = scenario$x_recent,
    xB = scenario$x_recent,
    dA = scenario$d_recent,
    dB = scenario$d_recent
  )

  old.parameters <- old.parameters[parameter.names]
  recent.parameters <- recent.parameters[parameter.names]

  geosse.parameters <- rbind(
    old = old.parameters,
    recent = recent.parameters
  )

  classe.parameters <- matrix(
    0,
    nrow = 2,
    ncol = length(diversitree:::default.argnames.classe(3))
  )

  colnames(classe.parameters) <- diversitree:::default.argnames.classe(3)
  rownames(classe.parameters) <- c("old", "recent")

  classe.parameters[1, ] <- divextra::pars.geosse2classe(
    geosse.parameters[1, ]
  )

  classe.parameters[2, ] <- divextra::pars.geosse2classe(
    geosse.parameters[2, ]
  )

  list(
    geosse = geosse.parameters,
    classe = classe.parameters
  )
}

simulate_tip_count <- function(scenario, total.time, seed, max.taxa) {
  simulation.parameters <- make_simulation_parameters(scenario)
  boundary.age <- total.time * scenario$boundary_fraction
  old.duration <- total.time - boundary.age

  set.seed(seed)
  simulated.root.state <- sample(1:3, size = 1)

  simulation.table <- divextra::make.tree.classe.td(
    pars.tb = simulation.parameters$classe,
    k = 3,
    max.taxa = max.taxa,
    max.t1 = old.duration,
    max.t2 = total.time,
    x0 = simulated.root.state,
    single.lineage = TRUE
  )

  if (is.null(simulation.table)) {
    return(0L)
  }

  simulated.time <- attr(simulation.table, "t.total")

  if (!isTRUE(all.equal(simulated.time, total.time, tolerance = 1e-8))) {
    return(as.integer(max.taxa))
  }

  tree <- tryCatch(
    divextra::table2tree(simulation.table),
    error = function(error) NULL
  )

  if (is.null(tree)) {
    return(0L)
  }

  as.integer(ape::Ntip(tree))
}

calibrate_scenario_duration <- function(
  scenario,
  calibration.replicates = 10,
  calibration.steps = 7,
  minimum.time = 1,
  initial.maximum.time = 10,
  maximum.time = 80
) {
  target.tips <- scenario$target_tips
  maximum.taxa <- max(2000, target.tips * 10)

  median_tips <- function(total.time, step.number) {
    tip.counts <- numeric()
    attempt.number <- 0L
    maximum.attempts <- calibration.replicates * 25

    while (
      length(tip.counts) < calibration.replicates &&
      attempt.number < maximum.attempts
    ) {
      attempt.number <- attempt.number + 1L
      seed <- 7000000 +
        scenario$scenario_number * 10000 +
        step.number * 100 +
        attempt.number

      tip.count <- simulate_tip_count(
        scenario = scenario,
        total.time = total.time,
        seed = seed,
        max.taxa = maximum.taxa
      )

      if (tip.count > 0) {
        tip.counts <- c(tip.counts, tip.count)
      }
    }

    if (length(tip.counts) < calibration.replicates) {
      stop(
        "Too few surviving calibration trees for scenario ",
        scenario$scenario_id,
        " at total.time = ", total.time
      )
    }

    median(tip.counts)
  }

  lower.time <- minimum.time
  upper.time <- initial.maximum.time
  upper.median <- median_tips(upper.time, 1)
  expansion.number <- 1

  while (upper.median < target.tips && upper.time < maximum.time) {
    lower.time <- upper.time
    upper.time <- min(upper.time * 1.5, maximum.time)
    expansion.number <- expansion.number + 1
    upper.median <- median_tips(upper.time, expansion.number)
  }

  if (upper.median < target.tips) {
    stop(
      "Could not reach the target tip count for scenario ",
      scenario$scenario_id,
      " before maximum.time = ", maximum.time
    )
  }

  final.median <- upper.median

  for (step.number in seq_len(calibration.steps)) {
    middle.time <- (lower.time + upper.time) / 2
    final.median <- median_tips(
      middle.time,
      expansion.number + step.number
    )

    if (final.median < target.tips) {
      lower.time <- middle.time
    } else {
      upper.time <- middle.time
    }
  }

  calibrated.time <- (lower.time + upper.time) / 2
  validation.median <- median_tips(
    calibrated.time,
    expansion.number + calibration.steps + 1
  )

  data.frame(
    total_time = calibrated.time,
    calibration_median_tips = validation.median,
    calibration_replicates = calibration.replicates,
    stringsAsFactors = FALSE
  )
}

check_benchmark_simulation <- function(
  simulation.table,
  tree,
  target.tips,
  epoch1.duration,
  tolerance = 1e-8
) {
  if (is.null(simulation.table)) {
    return(list(accepted = FALSE, reason = "total_extinction"))
  }

  if (is.null(tree)) {
    return(list(accepted = FALSE, reason = "tree_conversion_failed"))
  }

  epoch1.time <- attr(simulation.table, "t.epoch1")

  if (
    is.null(epoch1.time) ||
    !is.finite(epoch1.time) ||
    epoch1.time < epoch1.duration - tolerance
  ) {
    return(list(accepted = FALSE, reason = "target_reached_in_epoch1"))
  }

  observed.tips <- ape::Ntip(tree)

  if (observed.tips != target.tips) {
    return(list(accepted = FALSE, reason = "wrong_number_of_tips"))
  }

  simulated.time <- attr(simulation.table, "t.total")

  if (
    is.null(simulated.time) ||
    !is.finite(simulated.time) ||
    simulated.time <= epoch1.duration
  ) {
    return(list(accepted = FALSE, reason = "invalid_total_time"))
  }

  tree.height <- max(ape::node.depth.edgelength(tree))
  boundary.age <- simulated.time - epoch1.duration

  if (
    !is.finite(tree.height) ||
    boundary.age <= 0 ||
    boundary.age >= tree.height
  ) {
    return(list(accepted = FALSE, reason = "boundary_outside_reconstructed_tree"))
  }

  list(
    accepted = TRUE,
    reason = "accepted",
    observed.tips = observed.tips,
    epoch1.time = epoch1.time,
    simulated.time = simulated.time,
    tree.height = tree.height,
    boundary.age = boundary.age,
    realized.boundary.fraction = boundary.age / tree.height
  )
}

simulate_benchmark_dataset <- function(
  scenario,
  replicate.number,
  maximum.attempts = 500
) {
  simulation.parameters <- make_simulation_parameters(scenario)
  reference.total.time <- scenario$total_time
  reference.boundary.age <-
    reference.total.time * scenario$boundary_fraction
  old.duration <- reference.total.time - reference.boundary.age
  target.tips <- as.integer(scenario$target_tips)

  seed <- scenario$scenario_number * 100000 + replicate.number
  set.seed(seed)
  rejection.reasons <- character()

  for (attempt.number in seq_len(maximum.attempts)) {
    simulated.root.state <- sample(1:3, size = 1)

    simulation.table <- divextra::make.tree.classe.td(
      pars.tb = simulation.parameters$classe,
      k = 3,
      max.taxa = target.tips,
      max.t1 = old.duration,
      max.t2 = Inf,
      x0 = simulated.root.state,
      single.lineage = TRUE
    )

    if (is.null(simulation.table)) {
      rejection.reasons <- c(rejection.reasons, "total_extinction")
      next
    }

    tree <- tryCatch(
      divextra::table2tree(simulation.table),
      error = function(error) NULL
    )

    checks <- check_benchmark_simulation(
      simulation.table = simulation.table,
      tree = tree,
      target.tips = target.tips,
      epoch1.duration = old.duration
    )

    if (!checks$accepted) {
      rejection.reasons <- c(rejection.reasons, checks$reason)
      next
    }

    reconstructed.root.state <- as.integer(tree$node.state[1])

    if (
      is.na(reconstructed.root.state) ||
      reconstructed.root.state < 1 ||
      reconstructed.root.state > 3
    ) {
      rejection.reasons <- c(rejection.reasons, "invalid_root_state")
      next
    }

    rejection.counts <- table(factor(
      rejection.reasons,
      levels = c(
        "total_extinction",
        "tree_conversion_failed",
        "target_reached_in_epoch1",
        "wrong_number_of_tips",
        "invalid_total_time",
        "boundary_outside_reconstructed_tree",
        "invalid_root_state"
      )
    ))

    return(list(
      seed = seed,
      attempt = attempt.number,
      table = simulation.table,
      tree = tree,
      simulated.root.state = simulated.root.state,
      reconstructed.root.state = reconstructed.root.state,
      epoch1.duration = old.duration,
      simulated.time = checks$simulated.time,
      boundary.age = checks$boundary.age,
      tree.height = checks$tree.height,
      realized.boundary.fraction = checks$realized.boundary.fraction,
      rejection.counts = rejection.counts,
      parameters = simulation.parameters
    ))
  }

  stop(
    "Could not simulate an acceptable dataset after ",
    maximum.attempts,
    " attempts for scenario ", scenario$scenario_id,
    ", replicate ", replicate.number
  )
}

fit_benchmark_model <- function(
  tree,
  boundary.age,
  model.file,
  model.name,
  eps.values = c(0.1, 0.4, 0.8)
) {
  model.configuration <- divextra::read_yaml_pars_td(model.file)
  model.configuration$epoch.times <- boundary.age

  likelihood <- divextra::make.classe.td(
    tree = tree,
    states = tree$tip.state,
    k = 3,
    n.epoch = 2,
    sampling.f = c(1, 1, 1),
    strict = FALSE,
    control = list(backend = "gslode")
  )

  constraints <- divextra::make_constraints_sse_td(model.configuration)
  constrained.likelihood <- diversitree::constrain(
    likelihood,
    formulae = constraints
  )

  fit.results <- vector("list", length(eps.values))

  for (start.number in seq_along(eps.values)) {
    eps <- eps.values[start.number]

    starting.point <- divextra::init.pars.classe_td(
      lik.const.td = constrained.likelihood,
      phy = tree,
      k = 3,
      n.epoch = 2,
      eps = eps
    )

    start.time <- proc.time()["elapsed"]

    fit <- tryCatch(
      diversitree::find.mle(
        constrained.likelihood,
        starting.point,
        condition.surv = TRUE,
        root = diversitree::ROOT.OBS,
        keep.func = FALSE,
        control = list(maxit = 20000)
      ),
      error = function(error) error
    )

    elapsed <- unname(proc.time()["elapsed"] - start.time)

    if (inherits(fit, "error")) {
      fit.results[[start.number]] <- list(
        model = model.name,
        start = start.number,
        eps = eps,
        success = FALSE,
        error = conditionMessage(fit),
        elapsed_seconds = elapsed
      )
      next
    }

    parameter.count <- length(fit$par)
    sample.size <- ape::Ntip(tree)
    aic <- 2 * parameter.count - 2 * fit$lnLik

    if (sample.size > parameter.count + 1) {
      aicc <- aic +
        2 * parameter.count * (parameter.count + 1) /
        (sample.size - parameter.count - 1)
    } else {
      aicc <- NA_real_
    }

    fit.results[[start.number]] <- list(
      model = model.name,
      start = start.number,
      eps = eps,
      success = TRUE,
      error = NA_character_,
      elapsed_seconds = elapsed,
      lnLik = unname(fit$lnLik),
      AIC = unname(aic),
      AICc = unname(aicc),
      convergence = fit$convergence,
      message = fit$message,
      counts = fit$counts,
      par = fit$par,
      par.full = fit$par.full
    )
  }

  successful <- which(vapply(
    fit.results,
    function(result) isTRUE(result$success),
    logical(1)
  ))

  if (length(successful) == 0) {
    return(list(
      model = model.name,
      starts = fit.results,
      best.start = NA_integer_,
      best = NULL
    ))
  }

  log.likelihoods <- vapply(
    fit.results[successful],
    function(result) result$lnLik,
    numeric(1)
  )

  best.position <- which.max(log.likelihoods)
  best.start <- successful[best.position]

  list(
    model = model.name,
    starts = fit.results,
    best.start = best.start,
    best = fit.results[[best.start]]
  )
}

extract_epoch_estimates <- function(fit.result) {
  output <- c(
    s_recent = NA_real_,
    d_recent = NA_real_,
    x_recent = NA_real_,
    s_old = NA_real_,
    d_old = NA_real_,
    x_old = NA_real_
  )

  if (is.null(fit.result$best)) {
    return(output)
  }

  full.parameters <- fit.result$best$par.full

  output["s_recent"] <- full.parameters["lambda111.1"]
  output["d_recent"] <- full.parameters["q13.1"]
  output["x_recent"] <- full.parameters["mu1.1"]
  output["s_old"] <- full.parameters["lambda111.2"]
  output["d_old"] <- full.parameters["q13.2"]
  output["x_old"] <- full.parameters["mu1.2"]

  output
}

package_version_record <- function() {
  data.frame(
    package = c("R", "divextra", "diversitree", "ape"),
    version = c(
      paste(R.version$major, R.version$minor, sep = "."),
      as.character(utils::packageVersion("divextra")),
      as.character(utils::packageVersion("diversitree")),
      as.character(utils::packageVersion("ape"))
    ),
    stringsAsFactors = FALSE
  )
}
