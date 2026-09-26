# Descendant likelihoods for a complete-census reconstructed tree. Marginal
# internal-node ASRs are deliberately not multiplied into these messages.
.classe_td_descendant_cache <- function(backbone, schedule, extinction.cache,
                                         control = list()) {
  defaults <- list(rtol = extinction.cache$control$rtol,
                   atol = extinction.cache$control$atol,
                   max.step = extinction.cache$control$max.step,
                   interpolation.tol = extinction.cache$control$interpolation.tol,
                   max.refinements = 25L, max.points = extinction.cache$control$max.points)
  if (!is.list(control) || (length(control) &&
      (is.null(names(control)) || anyDuplicated(names(control)) ||
       any(!names(control) %in% names(defaults)))))
    stop("Unknown descendant-likelihood numerical control")
  control <- utils::modifyList(defaults, control, keep.null = TRUE)
  if (any(!vapply(control, function(z) is.numeric(z) && length(z) == 1L &&
                  is.finite(z) && z > 0, logical(1))))
    stop("Descendant controls must be finite positive numbers")
  k <- schedule$k
  tree <- backbone$tree
  ntip <- length(tree$tip.label)
  outgoing <- split(seq_len(nrow(tree$edge)), tree$edge[, 1])
  branches <- vector("list", nrow(tree$edge))
  scores <- vector("list", ntip + tree$Nnode)
  scales <- numeric(length(scores))
  rates <- extinction.cache$rates
  for (tip in seq_len(ntip)) {
    v <- backbone$tip.probs[, tip]
    scores[[tip]] <- v / sum(v)
    scales[tip] <- log(sum(v))
  }
  integrate.branch <- function(edge, initial, inherited) {
    lo <- backbone$edge.table$child.age[edge]
    hi <- backbone$edge.table$parent.age[edge]
    cuts <- sort(unique(c(lo, schedule$boundaries[
      schedule$boundaries > lo & schedule$boundaries < hi], hi)))
    times <- lo
    vals <- matrix(c(initial, 0), nrow = 1)
    current <- c(initial, 0)
    diagnostics <- list()
    if (hi > lo) for (segment in seq_len(length(cuts) - 1L)) {
      a <- cuts[segment]; b <- cuts[segment + 1L]
      epoch <- .classe_td_epoch((a + b) / 2, schedule$boundaries)
      rr <- rates[[epoch]]
      totals <- vapply(rr, function(r) sum(r$lambda, r$mu, r$q), numeric(1))
      linear <- diag(-totals,k)
      for (i in seq_len(k)) linear[i,rr[[i]]$q.to] <- rr[[i]]$q
      parents <- rep(seq_len(k),vapply(rr,function(r) length(r$lambda),integer(1)))
      pairs <- do.call(rbind,lapply(rr,`[[`,"lambda.ij"))
      birth.coefficient <- matrix(0,k,length(parents))
      if (length(parents)) birth.coefficient[cbind(parents,seq_along(parents))] <-
        unlist(lapply(rr,`[[`,"lambda"),use.names=FALSE)
      j <- pairs[,1]; l <- pairs[,2]
      rhs <- function(t, y, parms) {
        p <- y[seq_len(k)]
        ee <- extinction.cache$E(t)
        z <- as.vector(linear %*% p + birth.coefficient %*%
                         (ee[j]*p[l]+p[j]*ee[l]))
        cscale <- sum(z)
        list(c(z - cscale * p, cscale))
      }
      step <- min(control$max.step, 0.02 / max(totals, 1e-12))
      n <- max(1L, ceiling((b - a) / step))
      if (!is.finite(n) || n > control$max.points)
        stop("Initial descendant grid exceeds max.points; adjust numerical controls")
      near <- exp(seq(log(min(b-a, 1e-12)), log(min(b-a, step)), length.out=100L))
      base.grid <- sort(unique(c(seq(a,b,length.out=n+1L),a+near)))
      base.grid <- base.grid[base.grid <= b]
      accepted <- FALSE
      for (refine in seq_len(control$max.refinements)) {
        midpoint <- (head(base.grid,-1L)+tail(base.grid,-1L))/2
        grid <- sort(unique(c(base.grid,midpoint)))
        if (length(grid) > control$max.points)
          stop("Descendant interpolation exceeds max.points; adjust numerical controls")
        sol <- unname(deSolve::ode(current, grid, rhs, parms = NULL,
                                   method = "lsoda", rtol = control$rtol,
                                   atol = control$atol, tcrit = b,
                                   maxsteps = 100000L)[, -1, drop = FALSE])
        if (nrow(sol) != length(grid) || any(!is.finite(sol)) ||
            any(sol[, seq_len(k), drop = FALSE] < -10 * control$atol))
          stop("Descendant likelihood ODE failed; tighten numerical controls")
        mid <- match(midpoint,grid)
        left <- match(head(base.grid,-1L),grid)
        right <- match(tail(base.grid,-1L),grid)
        truth <- sol[mid,,drop=FALSE]
        interpolated <- (sol[left,,drop=FALSE]+sol[right,,drop=FALSE])/2
        errors <- abs(truth-interpolated)
        pp <- truth[,seq_len(k),drop=FALSE]
        pl <- interpolated[,seq_len(k),drop=FALSE]
        relative <- pp > 100*control$atol
        logerror <- matrix(0,nrow(pp),k)
        logerror[relative] <- abs(log(pp[relative])-log(pl[relative]))
        # Below the ODE's absolute-error scale, relative accuracy is not
        # attainable by refining interpolation. Use the solver's documented
        # absolute and relative tolerances in addition to interpolation error.
        logbudget <- control$interpolation.tol +
          10 * control$rtol + 10 * control$atol / pmax(pp, pl, control$atol)
        bad <- apply(errors,1,max)>control$interpolation.tol |
          apply(logerror > logbudget,1,any)
        if (!any(bad)) { accepted <- TRUE; break }
        base.grid <- sort(unique(c(base.grid,midpoint[bad])))
      }
      if (!accepted) stop(sprintf(paste0("Descendant interpolation did not converge on edge %d, ",
        "epoch %d: absolute error %.4g, log error %.4g; refine numerical controls"),
        edge,epoch,max(errors),max(logerror)))
      diagnostics[[segment]] <- data.frame(epoch=epoch, points=length(grid),
        refinements=refine-1L, max.absolute.error=max(errors),
        max.log.error=max(logerror),
        solver.limited.entries=sum(logerror > control$interpolation.tol))
      sol[, seq_len(k)][sol[, seq_len(k)] < 0] <- 0
      times <- c(times, grid[-1])
      vals <- rbind(vals, sol[-1, , drop = FALSE])
      current <- sol[nrow(sol), ]
    }
    list(times = times, probabilities = vals[, seq_len(k), drop = FALSE],
         log.scale = vals[, k + 1L] + inherited,
         diagnostics=if(length(diagnostics)) do.call(rbind,diagnostics) else NULL)
  }
  # Ascending node ages guarantee that daughter messages are ready first.
  nodes <- ntip + seq_len(tree$Nnode)
  for (node in nodes[order(backbone$age[nodes])]) {
    edges <- outgoing[[as.character(node)]]
    for (edge in edges) {
      child <- tree$edge[edge, 2]
      branches[[edge]] <- integrate.branch(edge, scores[[child]], scales[child])
    }
    v <- lapply(edges, function(e) tail(branches[[e]]$probabilities, 1L)[1, ])
    rr <- rates[[.classe_td_epoch(backbone$age[node], schedule$boundaries)]]
    d <- vapply(seq_len(k), function(i) {
      r <- rr[[i]]
      if (!length(r$lambda)) return(0)
      j <- r$lambda.ij[, 1]; l <- r$lambda.ij[, 2]
      sum(r$lambda * (v[[1]][j] * v[[2]][l] + v[[1]][l] * v[[2]][j]) / 2)
    }, numeric(1))
    if (!is.finite(sum(d)) || sum(d) <= 0)
      stop("Observed node has zero descendant likelihood under the supplied model")
    scores[[node]] <- d / sum(d)
    scales[node] <- log(sum(d)) + sum(vapply(edges, function(e)
      tail(branches[[e]]$log.scale, 1L), numeric(1)))
  }
  interpolators <- lapply(branches,function(b) {
    vv <- cbind(b$probabilities,b$log.scale)
    lapply(seq_len(k+1L),function(i) {
      if (length(b$times)==1L) {
        value <- vv[1,i]
        return(function(age) value)
      }
      stats::approxfun(b$times,vv[,i],ties="ordered",rule=1)
    })
  })
  interpolate <- function(age, edge) {
    b <- branches[[edge]]
    if (age < b$times[1] || age > tail(b$times, 1L))
      stop("Descendant likelihood requested outside its branch")
    vapply(interpolators[[edge]],function(f) f(age),numeric(1))
  }
  D <- function(age, edge) interpolate(age, edge)[seq_len(k)]
  logD <- function(age, state, edge) {
    ff <- interpolators[[edge]]
    log(ff[[state]](age)) + ff[[k+1L]](age)
  }
  list(root = scores[[backbone$root]], root.log.scale = scales[backbone$root],
       D = D, logD = logD, branches = branches, control = control,
       child.edges = function(node) outgoing[[as.character(node)]],
       child.scores = function(node) lapply(outgoing[[as.character(node)]],
         function(e) tail(branches[[e]]$probabilities, 1L)[1, ]))
}

.classe_td_draw_backbone_event <- function(age, state, stop.age, map, edge,
                                            schedule, Ecache, Dcache) {
  start <- Dcache$logD(age, state, edge)
  if (!is.finite(start)) stop("Cannot simulate from a zero-likelihood backbone state")
  rstart <- Ecache$integrated.rate(age, state)
  hazard <- function(u) rstart - Ecache$integrated.rate(u, state) +
    start - Dcache$logD(u, state, edge)
  target <- stats::rexp(1)
  end.hazard <- hazard(stop.age)
  if (is.na(end.hazard) || end.hazard < -1e-6)
    stop("Invalid backbone no-event probability; refine descendant controls")
  if (target >= end.hazard || age == stop.age)
    return(list(type = "end", age = stop.age, state = state,
                map = .append_map(map, state, age - stop.age, schedule$state.labels)))
  lower <- stop.age; upper <- age
  for (iteration in seq_len(200L)) {
    middle <- (lower + upper) / 2
    h <- hazard(middle)
    if (is.na(h) || h < -1e-6) stop("Invalid backbone conditional survival")
    if (h > target) lower <- middle else upper <- middle
    if (upper - lower <= Ecache$control$root.tol * max(upper, .Machine$double.eps)) break
  }
  event.age <- (lower + upper) / 2
  if (event.age <= stop.age || event.age >= age)
    stop("Backbone event age outside representable branch interval")
  epoch <- .classe_td_epoch(event.age, schedule$boundaries)
  r <- Ecache$rates[[epoch]][[state]]
  ee <- Ecache$E(event.age); dd <- Dcache$D(event.age, edge)
  pairs <- r$lambda.ij
  j <- pairs[, 1]; l <- pairs[, 2]
  weights <- c(r$lambda * dd[j] * ee[l], r$lambda * dd[l] * ee[j],
               r$q * dd[r$q.to])
  if (any(!is.finite(weights)) || sum(weights) <= 0)
    stop("Backbone conditional event weights are zero or invalid")
  z <- sample.int(length(weights), 1L, prob = weights)
  nb <- length(r$lambda)
  ans <- list(age = event.age, state = state, epoch = epoch,
              map = .append_map(map, state, age - event.age, schedule$state.labels))
  if (z <= 2L * nb) {
    ans$type <- "speciation"
    ans$pair <- if (z <= nb) pairs[z, ] else rev(pairs[z - nb, ])
  } else {
    ans$type <- "transition"
    ans$to <- r$q.to[z - 2L * nb]
  }
  ans
}
