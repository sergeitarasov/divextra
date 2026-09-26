#' Extract dated state segments from a ClaSSE stochastic map
#'
#' @param sim One result from [make.simmap.classe.td()] or
#'   [keep.lineage.simmap()].
#' @param tips Optional observed tip labels defining the focus. Applies
#'   [keep.lineage.simmap()] before extracting segments.
#' @return A data frame with lineage ID, parent and child graph IDs, whether
#'   the lineage is generated, state, and start/end ages before the present.
#'   Each row is one constant-state segment, not necessarily a whole species.
#' @details Uses exact numeric stochastic maps, not rounded history strings.
#'   Extinct tips retain their actual death ages. Ages use the tree's units.
#' @export
lineage.intervals.simmap <- function(sim, tips = NULL) {
  if (!is.null(tips)) sim <- keep.lineage.simmap(sim, tips)
  lin <- sim$lineages
  if (is.null(lin) || !nrow(lin) || is.null(sim$tree$maps) || is.null(sim$node.map))
    stop("sim must contain a mapped tree, lineages and node.map")
  ids <- setNames(sim$node.map$phylo.node, sim$node.map$graph.node)
  keys <- paste(sim$tree$edge[, 1], sim$tree$edge[, 2])
  idx <- match(paste(ids[lin$parent], ids[lin$child]), keys)
  rows <- lapply(seq_len(nrow(lin)), function(i) {
    x <- lin[i, ]
    if (!isTRUE(x$root.edge) && is.na(idx[i])) stop("Lineage missing from tree")
    m <- if (isTRUE(x$root.edge)) sim$tree$root.map else sim$tree$maps[[idx[i]]]
    if (!length(m) || is.null(names(m)) || any(!is.finite(m)) || any(m < 0) ||
        abs(sum(m) - (x$start.age - x$end.age)) > 1e-7)
      stop("Invalid mapped durations for lineage ", x$lineage.id)
    ends <- x$start.age - cumsum(m)
    ends[length(ends)] <- x$end.age
    starts <- c(x$start.age, head(ends, -1L))
    out <- data.frame(lineage.id = x$lineage.id, parent = x$parent, child = x$child,
      generated = x$generated, state = names(m), start.age = starts,
      end.age = ends, stringsAsFactors = FALSE)
    out[out$start.age > out$end.age, , drop = FALSE]
  })
  ans <- do.call(rbind, rows)
  rownames(ans) <- NULL
  ans
}

#' Focal-state lineage counts through time across stochastic maps
#'
#' @param sims One simulation, a list of simulations from
#'   [make.simmap.classe.td()], or a list of focused results.
#' @param times Nonnegative finite ages before present, in tree units (e.g. Ma).
#' @param states Exact state labels to pool, e.g. `c("S", "A.S", "S.R")`.
#'   Each lineage counts once if its current state belongs to this set.
#' @param tips Optional observed tip labels. The same focus selection is
#'   applied independently to every simulation using [keep.lineage.simmap()].
#' @param probs Two increasing quantile probabilities, default 0.025 and 0.975.
#' @param quantile.type Quantile algorithm passed to [stats::quantile()], an
#'   integer from 1 to 9. Default 1 selects observed counts; 7 interpolates.
#' @return A `classe_td_ltt` list containing `counts` (time by simulation),
#'   `summary` (time, mean, median, lower, upper, n.maps), `segments` (focal-state
#'   intervals with simulation IDs), `times`, `states`, `probs`, and `tips`.
#' @details
#' Counts living branches occupying any focal state, regardless of eventual
#' survival or extinction. At exact events uses the younger-side convention:
#' daughters count at birth, extinct lineages do not count at death, and state
#' transitions use the new state. Present-day counts use terminal tip states.
#' Before a map's root counts are zero. All maps, including zero-count maps,
#' contribute to every summary. Inputs must use a common age scale and state
#' coding. Quantiles are pointwise empirical map-distribution intervals
#' (type 1 by default), not confidence intervals for the mean or simultaneous bands.
#' With fixed parameters they do not include parameter uncertainty. For the
#' ASR-weighted mapper they describe that simulation distribution, not an exact
#' ClaSSE posterior. No extra conditioning or rejection is performed here.
#' @export
lineages.through.time.simmap <- function(sims, times, states, tips = NULL,
                                         probs = c(.025, .975), quantile.type = 1L) {
  .classe_td_check_quantile_type(quantile.type)
  if (!is.null(sims$tree)) sims <- list(sims)
  if (!is.list(sims) || !length(sims)) stop("sims must contain at least one map")
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) || any(times < 0))
    stop("times must be nonnegative finite ages")
  if (!is.character(states) || !length(states) || anyNA(states) || anyDuplicated(states))
    stop("states must be unique state labels")
  if (!is.numeric(probs) || length(probs) != 2L || any(!is.finite(probs)) ||
      probs[1] < 0 || probs[2] > 1 || probs[1] >= probs[2])
    stop("probs must contain two increasing probabilities in [0, 1]")
  times <- sort(unique(times))
  counts <- matrix(0L, length(times), length(sims),
                   dimnames = list(as.character(times), paste0("sim", seq_along(sims))))
  segments <- vector("list", length(sims))
  for (i in seq_along(sims)) {
    sim <- sims[[i]]
    if (!all(states %in% colnames(sim$tree$mapped.edge)))
      stop("Unknown focal state in simulation ", i)
    if (!is.null(tips)) sim <- keep.lineage.simmap(sim, tips)
    z <- lineage.intervals.simmap(sim)
    z <- z[z$state %in% states, , drop = FALSE]
    # Number started minus number ended, moving from present into the past.
    # left.open=TRUE gives strict '< time', hence end < time <= start.
    counts[, i] <- findInterval(times, sort(z$end.age), left.open = TRUE) -
      findInterval(times, sort(z$start.age), left.open = TRUE)
    if (any(times == 0)) counts[times == 0, i] <-
      sum(sim$tips$age == 0 & sim$tips$state.label %in% states)
    z$simulation <- rep.int(i, nrow(z))
    segments[[i]] <- z
  }
  summary <- data.frame(time = times, mean = rowMeans(counts),
    median = apply(counts, 1, stats::median),
    lower = apply(counts, 1, stats::quantile, probs = probs[1], type = quantile.type),
    upper = apply(counts, 1, stats::quantile, probs = probs[2], type = quantile.type),
    n.maps = length(sims), row.names = NULL)
  structure(list(counts = counts, summary = summary,
    segments = do.call(rbind, segments), times = times, states = states,
    probs = probs, tips = tips, quantile.type = quantile.type), class = "classe_td_ltt")
}

#' Plot a focal-state lineage-through-time distribution
#'
#' @param x Output of [lineages.through.time.simmap()].
#' @param type `"interval"` for median and quantile band, `"maps"` for
#'   translucent count curves for individual maps, or `"segments"` for
#'   individual focal-state branch segments in one map.
#' @param simulation Map index for `type = "segments"`.
#' @param col Line/segment colour.
#' @param alpha Transparency for the band or individual-map curves.
#' @param xlab,ylab Axis labels. Time is shown oldest to youngest.
#' @param center Central curve: `"median"` (default) or `"mean"`.
#' @param interval Band: `"auto"` uses map quantiles for the median and a
#'   confidence interval for the mean; `"distribution"` uses map quantiles;
#'   `"mean-ci"` requires `center="mean"`. Only used for `type="interval"`.
#' @param level Confidence level for the mean interval, default 0.95.
#' @param quantile.type Optional quantile algorithm (1 to 9) overriding the
#'   distribution band using stored counts, without changing `x`. NULL retains
#'   the stored bounds. Ignored for mean confidence intervals and other plot types.
#' @param ... Further arguments to the initial base plot, e.g. `main`.
#' @return Invisibly returns `x`.
#' @details Lines between sampled ages are visual interpolation, not additional
#'   counts. Segment plots use exact event ages; each row is a branch lineage,
#'   not a continuous species identity across speciation events.
#' Mean intervals are pointwise Student-t intervals using sd/sqrt(n), truncated
#' below at zero. At least two independent maps are required. These approximate
#' intervals measure Monte Carlo uncertainty in the mean, not variation among
#' histories or parameter uncertainty; small skewed samples can give poor coverage.
#' Pooled maps are treated as independent draws from one mixture, not as a
#' stratified sample with fixed model allocations. Distribution bands retain
#' `x$probs`; `level` affects only mean confidence intervals. Segment plots are
#' unaffected by these options.
#' @export
plot.classe_td_ltt <- function(x, type = c("interval", "maps", "segments"),
                               simulation = 1L, col = "#0072B2", alpha = .15,
                               xlab = "Age before present", ylab = NULL,
                               center = c("median", "mean"),
                               interval = c("auto", "distribution", "mean-ci"),
                               level = .95, quantile.type = NULL, ...) {
  if (!is.null(quantile.type)) .classe_td_check_quantile_type(quantile.type)
  type <- match.arg(type)
  center <- match.arg(center)
  interval <- match.arg(interval)
  if (interval == "auto") interval <- if (center == "mean") "mean-ci" else "distribution"
  s <- x$summary
  if (type == "segments") {
    if (length(simulation) != 1L || is.na(simulation) ||
        !simulation %in% seq_len(ncol(x$counts))) stop("Invalid simulation index")
    z <- x$segments[x$segments$simulation == simulation, , drop = FALSE]
    ids <- unique(z$lineage.id)
    if (is.null(ylab)) ylab <- "Branch lineage"
    graphics::plot(NA, xlim = rev(range(c(x$times, z$start.age, z$end.age))),
      ylim = c(.5, max(1, length(ids)) + .5), xlab = xlab, ylab = ylab, yaxt = "n", ...)
    if (nrow(z)) {
      y <- match(z$lineage.id, ids)
      graphics::segments(z$start.age, y, z$end.age, y, col = col, lwd = 3)
      graphics::axis(2, at = seq_along(ids), labels = ids, las = 1)
    }
  } else {
    central <- s[[center]]
    lower <- s$lower; upper <- s$upper
    if (type == "interval" && interval == "distribution" && !is.null(quantile.type)) {
      lower <- apply(x$counts,1,stats::quantile,probs=x$probs[1],type=quantile.type)
      upper <- apply(x$counts,1,stats::quantile,probs=x$probs[2],type=quantile.type)
    }
    if (type == "interval" && interval == "mean-ci") {
      if (center != "mean") stop("mean-ci requires center='mean'")
      ci <- .classe_td_ltt_mean_ci(x$counts, level)
      lower <- ci[,1]; upper <- ci[,2]
    }
    if (is.null(ylab)) ylab <- "Number of focal-state lineages"
    top <- if (type == "maps") max(x$counts) else max(upper, central)
    graphics::plot(s$time, central, type = "n", xlim = rev(range(s$time)),
      ylim = c(0, max(1, top)), xlab = xlab, ylab = ylab, ...)
    if (type == "interval") {
      graphics::polygon(c(s$time, rev(s$time)), c(lower, rev(upper)),
                        col = grDevices::adjustcolor(col, alpha.f = alpha), border = NA)
      graphics::segments(s$time, lower, s$time, upper,
                         col = grDevices::adjustcolor(col, alpha.f = alpha))
    } else {
      graphics::matlines(x$times, x$counts, col = grDevices::adjustcolor(col, alpha.f = alpha),
                        lty = 1, lwd = 1)
    }
    graphics::lines(s$time, central, col = col, lwd = 2)
    graphics::points(s$time, central, col = col, pch = 16, cex = .4)
  }
  invisible(x)
}

.classe_td_ltt_mean_ci <- function(counts, level) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must be between 0 and 1")
  n <- ncol(counts)
  if (n < 2L) stop("Mean confidence intervals require at least two independent maps")
  average <- rowMeans(counts)
  margin <- stats::qt((1+level)/2, df=n-1) * apply(counts, 1, stats::sd) / sqrt(n)
  cbind(lower=pmax(0, average-margin), upper=average+margin)
}

.classe_td_check_quantile_type <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      !x %in% 1:9) stop("quantile.type must be an integer from 1 to 9")
  invisible(x)
}
