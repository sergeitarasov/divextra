# Shared extinction probabilities for the unrestricted branching process.
# E means biological extinction, not absence from a sampled tree: E(0) = 0.
# The interpolation grid is adaptively checked against new ODE evaluations.
# These checks are numerical convergence diagnostics, not rigorous error bounds.
.classe_td_extinction_cache <- function(schedule, max.age, control = list()) {
  defaults <- list(rtol = 1e-10, atol = 1e-13, max.step = 0.01,
                   root.tol = 1e-10, interpolation.tol = 1e-6,
                   max.points = 250000)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyDuplicated(names(control)) ||
       any(!names(control) %in% names(defaults)))))
    stop(paste("Unknown extinction control; use rtol, atol, max.step,",
               "root.tol, interpolation.tol, max.points"))
  control <- utils::modifyList(defaults, control, keep.null = TRUE)
  if (any(!vapply(control, function(x) is.numeric(x) && length(x) == 1L &&
                  is.finite(x) && x > 0, logical(1))))
    stop("Extinction controls must be finite positive numbers")
  if (length(max.age) != 1L || !is.finite(max.age) || max.age <= 0)
    stop("max.age must be positive and finite")
  k <- schedule$k
  rates <- lapply(schedule$arrays, function(a)
    lapply(seq_len(k), function(i) .classe_td_event_rates(a, i)))
  totals <- vapply(rates, function(rr)
    vapply(rr, function(r) sum(r$lambda, r$mu, r$q), numeric(1)), numeric(k))
  totals <- matrix(totals, nrow = k)
  boundaries <- schedule$boundaries
  integral <- function(age, state) {
    lower <- c(0, boundaries)
    upper <- c(boundaries, Inf)
    sum(totals[state, ] * pmax(0, pmin(age, upper) - lower))
  }
  cuts <- c(0, boundaries[boundaries > 0 & boundaries < max.age], max.age)
  values <- matrix(0, 1L, k)
  times <- 0
  initial <- numeric(k)
  reachable <- rep(FALSE, k)
  support <- list()
  diagnostics <- list()
  for (segment in seq_len(length(cuts) - 1L)) {
    lo <- cuts[segment]
    hi <- cuts[segment + 1L]
    epoch <- .classe_td_epoch((lo + hi) / 2, boundaries)
    rr <- rates[[epoch]]
    # Structural support: a lineage can die directly, change to a supported
    # state, or split into two supported states. Carry support across epochs.
    repeat {
      previous <- reachable
      for (i in seq_len(k)) {
        r <- rr[[i]]
        births <- if (length(r$lambda))
          any(reachable[r$lambda.ij[, 1]] & reachable[r$lambda.ij[, 2]]) else FALSE
        reachable[i] <- reachable[i] || r$mu > 0 ||
          any(reachable[r$q.to]) || births
      }
      if (identical(previous, reachable)) break
    }
    support[[segment]] <- reachable
    rhs <- function(t, y, parms) {
      dy <- vapply(seq_len(k), function(i) {
        r <- rr[[i]]
        births <- if (length(r$lambda))
          sum(r$lambda * y[r$lambda.ij[, 1]] * y[r$lambda.ij[, 2]]) else 0
        r$mu - totals[i, epoch] * y[i] + sum(r$q * y[r$q.to]) + births
      }, numeric(1))
      dy[!reachable] <- 0
      list(dy)
    }
    span <- hi - lo
    # Geometric points resolve E near zero, including states whose first
    # possible route to extinction becomes available at an epoch boundary.
    near <- exp(seq(log(min(span, 1e-12)),
                    log(min(span, control$max.step)), length.out = 160L))
    if (ceiling(span / control$max.step) + length(times) > control$max.points)
      stop("Extinction grid exceeds max.points; increase max.points or max.step")
    # Resolve fast initial transients without a globally microscopic time step.
    rate <- max(totals[, epoch])
    transient <- if (rate > 0)
      expm1(seq(0, log1p(min(rate * span, 20)), length.out = 160L)) / rate else 0
    grid <- sort(unique(c(lo, hi, lo + near, lo + transient,
                         seq(lo, hi, length.out = ceiling(span / control$max.step) + 1L))))
    grid <- grid[grid >= lo & grid <= hi]
    solve.grid <- function(at) {
      sol <- deSolve::ode(initial, at, rhs, parms = NULL, method = "lsoda",
                          rtol = control$rtol, atol = control$atol,
                          maxsteps = 100000L)
      ans <- unname(sol[, -1, drop = FALSE])
      if (nrow(ans) != length(at) || any(!is.finite(ans)) ||
          any(ans < -10 * control$atol) ||
          any(ans > 1 + 10 * (control$atol + control$rtol)))
        stop("Extinction ODE failed its probability checks; tighten solver controls")
      # Only clip integration roundoff outside [0,1], never floor small E.
      ans[ans < 0] <- 0
      ans[ans > 1] <- 1
      ans[, !reachable] <- 0
      ans
    }
    for (refinement in seq_len(30L)) {
      n <- length(grid)
      if (n + length(times) - 1L > control$max.points)
        stop("Adaptive extinction grid exceeds max.points; refine solver controls or increase max.points")
      left <- grid[-n]
      width <- diff(grid)
      fraction <- rep(c(.25, .5, .75), each = n - 1L)
      interval <- rep(seq_len(n - 1L), 3L)
      probe <- left[interval] + width[interval] * fraction
      all.times <- sort(unique(c(grid, probe)))
      all.values <- solve.grid(all.times)
      yy <- all.values[match(grid, all.times), , drop = FALSE]
      truth <- all.values[match(probe, all.times), , drop = FALSE]
      linear <- yy[interval, , drop = FALSE] * (1 - fraction) +
        yy[interval + 1L, , drop = FALSE] * fraction
      absolute.error <- abs(truth - linear)
      # Relative precision cannot be certified below the ODE's absolute
      # resolution. Expose that threshold rather than silently flooring E.
      resolved <- truth > 100 * control$atol
      log.error <- matrix(0, nrow(truth), k)
      log.error[resolved] <- abs(log(truth[resolved]) - log(linear[resolved]))
      bad <- apply(absolute.error > control$interpolation.tol |
                     log.error > control$interpolation.tol, 1L, any)
      if (!any(bad)) break
      additions <- probe[interval %in% unique(interval[bad])]
      updated <- sort(unique(c(grid, additions)))
      if (identical(updated, grid) || refinement == 30L)
        stop("Extinction interpolation failed to converge; tighten ODE controls")
      grid <- updated
    }
    diagnostics[[segment]] <- list(epoch = epoch, points = length(grid),
      refinements = refinement - 1L, max.absolute.error = max(absolute.error),
      max.log.error = max(log.error), relative.check.above = 100 * control$atol)
    initial <- yy[nrow(yy), ]
    times <- c(times, grid[-1])
    values <- rbind(values, yy[-1, , drop = FALSE])
  }
  interpolators <- lapply(seq_len(k), function(i)
    stats::approxfun(times, values[, i], rule = 1, ties = "ordered"))
  E <- function(age) {
    if (length(age) != 1L || !is.finite(age) || age < 0 || age > max.age)
      stop("Extinction probability requested outside the cached age range")
    vapply(interpolators, function(f) f(age), numeric(1))
  }
  Ei <- function(age, state) interpolators[[state]](age)
  possible <- function(age, state) {
    segment <- sum(age > cuts[-length(cuts)])
    segment > 0L && support[[segment]][state]
  }
  list(E = E, Ei = Ei, integrated.rate = integral, times = times, values = values,
       possible = possible, rates = rates, control = control, max.age = max.age,
       diagnostics = diagnostics)
}

# Invert conditional no-event survival rather than integrating mu/E near zero.
# This avoids singular conditional rates at the extinction deadline.
.classe_td_draw_extinct_event <- function(age, state, map, schedule, cache) {
  e.start <- cache$Ei(age, state)
  if (!cache$possible(age, state))
    stop("Cannot condition side lineage on extinction: extinction is structurally impossible")
  if (!is.finite(e.start) || e.start <= 0)
    stop("Extinction probability is below numerical resolution; tighten solver controls")
  r.start <- cache$integrated.rate(age, state)
  target <- stats::rexp(1)
  hazard <- function(u) {
    eu <- cache$Ei(u, state)
    if (eu <= 0) return(Inf)
    r.start - cache$integrated.rate(u, state) + log(e.start) - log(eu)
  }
  # The zero-age endpoint has infinite hazard. Bisect directly so no infinite
  # endpoint is passed to uniroot and keep the root strictly above the deadline.
  lower <- 0
  upper <- age
  for (iteration in seq_len(200L)) {
    middle <- (lower + upper) / 2
    if (middle == lower || middle == upper) break
    h <- hazard(middle)
    if (is.na(h) || h < -1e-7)
      stop("Invalid conditional survival; refine extinction interpolation")
    if (h > target) lower <- middle else upper <- middle
    # Relative tolerance also resolves extremely short histories near time zero.
    if (upper - lower <= cache$control$root.tol * max(upper, .Machine$double.eps)) break
  }
  event.age <- (lower + upper) / 2
  if (event.age <= 0 || event.age >= age)
    stop("Conditional event time is outside representable open age interval")
  epoch <- .classe_td_epoch(event.age, schedule$boundaries)
  r <- cache$rates[[epoch]][[state]]
  e <- cache$E(event.age)
  birth.weights <- if (length(r$lambda))
    r$lambda * e[r$lambda.ij[, 1]] * e[r$lambda.ij[, 2]] else numeric()
  weights <- c(birth.weights, r$mu, r$q * e[r$q.to])
  if (any(!is.finite(weights)) || sum(weights) <= 0)
    stop("Conditional event weights are numerically zero or invalid")
  # Preserve epoch crossings in chronological map order (adjacent equal state
  # segments are merged by the existing map helper).
  map <- .append_map(map, state, age - event.age, schedule$state.labels)
  z <- sample.int(length(weights), 1L, prob = weights)
  nl <- length(birth.weights)
  ans <- list(age = event.age, epoch = epoch, state = state, map = map)
  if (z <= nl) {
    ans$type <- "speciation"
    ans$pair <- r$lambda.ij[z, ]
  } else if (z == nl + 1L) {
    ans$type <- "extinction"
  } else {
    ans$type <- "transition"
    ans$to <- r$q.to[z - nl - 1L]
  }
  ans
}
