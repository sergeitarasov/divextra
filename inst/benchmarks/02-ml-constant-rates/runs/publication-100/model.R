# r-env may override R_LIBS_USER; explicitly prepend the project library.
project.library <- "/projappl/project_2004996/geosse-benchmark/Rlib/4.6.1"
if (dir.exists(project.library)) .libPaths(c(project.library, .libPaths()))
library(divextra)

# Identical likelihood settings on the local computer and Roihu.
make.likelihood <- function(z, model) {
    cfg <- read_yaml_pars_td("model-episodic.yml")
    cfg$epoch.times <- z$recent.time
    raw <- make.classe.td(z$tree, z$tree$tip.state, k = 3, n.epoch = 2,
      strict = FALSE, sampling.f = c(1, 1, 1), control = list(tol = 1e-10, eps = 1e-8))
    raw <- diversitree::constrain(raw, formulae = make_constraints_sse_td(cfg))
    stopifnot(identical(diversitree::argnames(raw),
      c("lambda111.1", "mu1.1", "q13.1", "lambda111.2", "mu1.2", "q13.2")))
    evaluate <- function(p) {
      full <- setNames(numeric(6), names(z$truth))
      for (rate in c("s", "d", "x")) {
        full[paste0(rate, "_recent")] <- p[paste0(rate, "_recent")]
        old.name <- if (z$links[rate]) paste0(rate, "_recent") else paste0(rate, "_old")
        full[paste0(rate, "_old")] <- p[old.name]
      }
      raw(full[c(1, 3, 2, 4, 6, 5)], root = diversitree::ROOT.OBS, condition.surv = TRUE)
    }
  evaluate
}

fit.one <- function(z, job) {
  evaluate <- make.likelihood(z, job$model)
  objective <- function(logp) {
    value <- tryCatch(suppressWarnings(evaluate(exp(logp))), error = function(e) -Inf)
    if (!is.finite(value)) return(1e50)
    -value
  }
  begin <- proc.time()["elapsed"]
  fit <- nlminb(log(job$initial), objective, lower = log(1e-10), upper = log(100),
    control = list(iter.max = 1000, eval.max = 3000, rel.tol = 1e-9))
  # Preserve the primary fit; separately test its near-zero-rate boundary.
  boundary <- NULL
  pars <- exp(fit$par)
  zero <- which(pars < 1e-5)
  free <- setdiff(seq_along(pars), zero)
  if (length(zero) > 0 && length(free) > 0) {
    boundary.objective <- function(logp) {
      candidate <- pars
      candidate[zero] <- 0
      candidate[free] <- exp(logp)
      value <- tryCatch(suppressWarnings(evaluate(candidate)), error = function(e) -Inf)
      if (!is.finite(value)) return(1e50)
      -value
    }
    bf <- tryCatch(nlminb(log(pars[free]), boundary.objective, lower = log(1e-10), upper = log(100),
      control = list(iter.max = 1000, eval.max = 3000, rel.tol = 1e-12)),
      error = function(e) list(par = log(pars[free]), objective = NA_real_,
        convergence = 99L, message = conditionMessage(e)))
    candidate <- pars
    candidate[zero] <- 0
    candidate[free] <- exp(bf$par)
    boundary <- list(par = candidate, lnLik = -bf$objective,
      convergence = bf$convergence, message = bf$message, fixed.zero = names(pars)[zero])
  }
  list(id = z$id, model = job$model, start = job$start,
    initial = job$initial, par = exp(fit$par), lnLik = -fit$objective,
    convergence = fit$convergence, message = fit$message,
    elapsed = unname(proc.time()["elapsed"] - begin), optimizer = fit, boundary = boundary)
}
