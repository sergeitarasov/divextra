#' Simulate ClaSSE with a stopping rule for each epoch
#'
#' @param pars.tb Numeric matrix, one row per epoch, oldest first. Columns use
#'   the standard diversitree ClaSSE order (or its complete parameter names).
#' @param k Number of states, at least two.
#' @param ntips Target number of living lineages per epoch. A scalar is recycled;
#'   Inf disables the target. Targets are exact first-hitting counts.
#' @param duration Maximum duration of each epoch, measured from its start,
#'   not a cumulative age. A scalar is recycled; Inf disables the time limit.
#' @param x0 Initial state, in 1:k.
#' @param max.taxa Global safety cap on living lineages, not an epoch target.
#' @param max.events Global safety cap on simulated events.
#' @return A lineage data frame compatible with table2simmap. Attributes include
#'   hist, k, sim.pars, t.total, info.epochs (boundary snapshots), epochs
#'   (realized forward times, ages before simulation end, counts and reasons),
#'   completed and stop.reason. Early termination is returned, never retried.
#' @details
#' Based on the Gillespie construction in make.tree.classe.td and diversitree.
#' The original function is unchanged. Any number of epochs is supported.
#' The first of the target count and duration limits ends an epoch. A target
#' already equal to the current count causes a zero-duration epoch; a lower
#' target requires deaths to reach it. Counts include all living lineages,
#' including those which may go extinct later, not reconstructed crossings.
#'
#' Tip stopping occurs immediately after the event that reaches the target,
#' without waiting for another birth. Terminal branches may therefore have zero
#' length. Extinction, an unreachable target with zero total event rate, or a
#' safety limit stops the whole simulation. At least one epoch limit must be
#' finite. Durations and ntips must have length one or the number of epochs.
#' Tip-dependent boundaries condition the sampling scheme; they do not guarantee
#' unbiased recovery under a fixed-boundary likelihood. In inference reverse the
#' parameter rows and use the realized ages, not requested durations.
#' Legacy table2tree assumes two epochs and is not a general converter for this
#' output; use table2simmap for mapped trees.
#' @export
#' @examples
#' p <- pars.geosse2classe(c(.1,.1,.1,.05,.05,.1,.1))
#' set.seed(10)
#' tb <- make.tree.classe.episodic(rbind(p,p), k=3, x0=1,
#'   ntips=c(100,500), duration=c(50,20))
#' attr(tb, "epochs")
make.tree.classe.episodic <- function(pars.tb, k, ntips=Inf, duration=Inf,
                                     x0=1, max.taxa=10000, max.events=1000000) {
  integer.scalar <- function(x, minimum) {
    is.numeric(x) && length(x)==1L && is.finite(x) &&
      x>=minimum && x==floor(x)
  }
  if (!integer.scalar(k,2)) stop("k must be an integer >= 2")
  if (!integer.scalar(x0,1) || x0>k) stop("x0 must be an integer in 1:k")
  if (!is.matrix(pars.tb) || !is.numeric(pars.tb) || nrow(pars.tb)<1L ||
      any(!is.finite(pars.tb)) || any(pars.tb<0))
    stop("pars.tb must be a finite nonnegative numeric matrix")
  expected <- diversitree:::default.argnames.classe(k)
  if (ncol(pars.tb)!=length(expected)) stop("Wrong number of ClaSSE parameters")
  if (!is.null(colnames(pars.tb))) {
    if (anyDuplicated(colnames(pars.tb)) || !setequal(colnames(pars.tb),expected))
      stop("Column names must be the complete ClaSSE parameter names")
    pars.tb <- pars.tb[,expected,drop=FALSE]
  }
  ne <- nrow(pars.tb)
  limits <- function(x, name, counts=FALSE) {
    if (!is.numeric(x) || !(length(x) %in% c(1L,ne)) ||
        anyNA(x) || any(x<if(counts) 1 else 0) ||
        (counts && any(is.finite(x) & x!=floor(x))))
      stop(name," must contain valid nonnegative durations or positive integer counts")
    rep(x,length.out=ne)
  }
  ntips <- limits(ntips,"ntips",TRUE)
  duration <- limits(duration,"duration")
  if (any(!is.finite(ntips) & !is.finite(duration)))
    stop("Each epoch needs a finite ntips or duration limit")
  if (!(integer.scalar(max.taxa,1) || identical(max.taxa,Inf)))
    stop("max.taxa must be a positive integer or Inf")
  if (!integer.scalar(max.events,1)) stop("max.events must be a positive integer")

  # Lineage and transition tables retain the layout of make.tree.classe.td.
  states <- x0
  parent <- start <- len <- 0
  extinct <- split <- FALSE
  alive <- 1L
  t <- 0
  events <- 0L
  hist <- list()
  snapshots <- list()
  epochs <- data.frame()
  destinations <- matrix(0L,k,k-1L)
  for(i in seq_len(k)) destinations[i,] <- setdiff(seq_len(k),i)
  completed <- FALSE
  for(epoch in seq_len(ne)) {
    begin <- t
    n.start <- length(alive)
    event.start <- events
    deadline <- begin+duration[epoch]
    rates <- diversitree:::inflate.pars.classe(pars.tb[epoch,],k)
    rates$lambda[is.na(rates$lambda)] <- 0
    q <- matrix(pars.tb[epoch,-seq_len(k*k*(k+1)/2+k)],k,k-1,byrow=TRUE)
    event.rates <- cbind(rowSums(rates$lambda),rates$mu,q)
    totals <- rowSums(event.rates)
    repeat {
      n <- length(alive)
      if(n==0) {reason <- "extinct"; break}
      if(n==ntips[epoch]) {reason <- "ntips"; break}
      if(t>=deadline) {reason <- "duration"; break}
      if(n>=max.taxa) {reason <- "max.taxa"; break}
      if(events>=max.events) {reason <- "max.events"; break}
      counts <- tabulate(states[alive],nbins=k)
      weights <- totals*counts
      total <- sum(weights)
      if(!is.finite(total)) stop("Total event rate overflow")
      if(total==0) {
        if(is.finite(deadline)) {
          len[alive] <- len[alive]+deadline-t
          t <- deadline
          reason <- "duration"
        } else reason <- "stalled"
        break
      }
      dt <- stats::rexp(1,total)
      if(t+dt>=deadline) {
        len[alive] <- len[alive]+deadline-t
        t <- deadline
        reason <- "duration"
        break
      }
      if(!is.finite(t+dt) || t+dt<=t) stop("Event time is not numerically representable")
      t <- t+dt
      len[alive] <- len[alive]+dt
      state <- sample.int(k,1,prob=weights)
      candidates <- alive[states[alive]==state]
      lineage <- candidates[sample.int(length(candidates),1)]
      type <- sample.int(k+1,1,prob=event.rates[state,])
      events <- events+1L
      if(type==1L) {
        j <- sample.int(k*k,1,prob=rates$lambda[state,,])
        daughters <- c((j-1L) %% k+1L,(j-1L) %/% k+1L)
        new <- length(states)+1:2
        split[lineage] <- TRUE
        states[new] <- daughters
        parent[new] <- lineage
        start[new] <- t
        len[new] <- 0
        split[new] <- extinct[new] <- FALSE
      } else if(type==2L) {
        extinct[lineage] <- TRUE
      } else {
        destination <- destinations[state,type-2L]
        states[lineage] <- destination
        hist[[length(hist)+1L]] <- c(lineage,t,state,destination)
      }
      alive <- which(!split & !extinct)
    }
    info <- data.frame(idx=seq_along(states),len=len,parent=parent,
      start=start,state=states,extinct=extinct,split=split)
    attr(info,"hist") <- .classe_td_history_table(hist,info)
    attr(info,"k") <- k
    snapshots[[epoch]] <- info
    epochs <- rbind(epochs,data.frame(epoch=epoch,start=begin,end=t,
      duration=t-begin,ntips.target=ntips[epoch],duration.limit=duration[epoch],
      n.start=n.start,n.end=length(alive),events=events-event.start,reason=reason))
    if(!reason %in% c("ntips","duration")) break
    completed <- epoch==ne
  }
  epochs$start.age <- t-epochs$start
  epochs$end.age <- t-epochs$end
  attr(info,"info.epochs") <- snapshots
  attr(info,"info.epoch1") <- snapshots[[1]]
  attr(info,"t.epoch1") <- epochs$end[1]
  attr(info,"t.total") <- t
  attr(info,"sim.pars") <- pars.tb
  attr(info,"epochs") <- epochs
  attr(info,"completed") <- completed
  attr(info,"stop.reason") <- reason
  info
}
