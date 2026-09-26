.classe_td_transition_events <- function(sim, event.types) {
  if (is.null(sim$tips$generated) || is.null(sim$lineages$generated) ||
      any(sim$tips$generated) || any(sim$lineages$generated))
    stop("Supply observed-only maps: call prune.unobserved.simmap first")
  z <- lineage.intervals.simmap(sim)
  empty <- data.frame(age=numeric(), from=character(), to=character(),
    event.type=character(), lineage.id=integer(), node=character())
  rows <- list(empty)
  groups <- split(z, z$lineage.id)
  for (g in groups) {
    if ("anagenetic" %in% event.types && nrow(g) > 1L) {
      j <- seq.int(2L, nrow(g))
      changed <- g$state[j] != g$state[j-1L]
      j <- j[changed]
      if (length(j)) rows[[length(rows)+1L]] <- data.frame(age=g$start.age[j],
        from=g$state[j-1L], to=g$state[j], event.type="anagenetic",
        lineage.id=g$lineage.id[j], node=g$parent[j])
    }
  }
  if ("cladogenetic" %in% event.types) {
    ev <- sim$events
    if (is.null(ev$from.label)) stop("Cladogenetic transitions require event history")
    ev <- ev[ev$event %in% c("observed_speciation", "hidden_speciation", "speciation"), ]
    for (g in groups) {
      j <- match(g$parent[1], ev$node)
      if (!is.na(j) && !is.na(ev$from.label[j]) && ev$from.label[j] != g$state[1])
        rows[[length(rows)+1L]] <- data.frame(age=g$start.age[1],
          from=ev$from.label[j], to=g$state[1], event.type="cladogenetic",
          lineage.id=g$lineage.id[1], node=g$parent[1])
    }
  }
  out <- do.call(rbind, rows)
  out <- out[order(out$age, decreasing=TRUE), , drop=FALSE]
  rownames(out) <- NULL
  out
}

.classe_td_probability <- function(successes, n, level) {
  if (!n) return(c(probability=NA_real_, lower=NA_real_, upper=NA_real_))
  ci <- stats::binom.test(successes, n, conf.level=level)$conf.int
  c(probability=successes/n, lower=ci[1], upper=ci[2])
}

#' Transition ages and probabilities on observed-only stochastic maps
#'
#' @param sims One result or a list of results from [prune.unobserved.simmap()].
#' @param breaks Increasing, nonnegative age-bin boundaries in tree units.
#'   Bins include the younger boundary and exclude the older, except that the
#'   last bin also includes its older boundary. Each event belongs to one bin.
#' @param to Exact destination state labels, e.g. `c("S", "A.S", "S.R")`.
#' @param from Exact source state labels. Default NULL selects every state
#'   outside `to`, so movement between destination states is not a new arrival.
#' @param event.types `"anagenetic"`, `"cladogenetic"`, or both (default).
#' @param level Interval coverage, default 0.95.
#' @return A `classe_td_transitions` list with:
#' \itemize{
#'   \item `events`: matching transitions, exact ages, type and map ID.
#'   \item `per.map`: total counts, counts by event type, oldest (first) and
#'     youngest (last) arrival ages, and `root.in.destination` (whether the
#'     root already occupied a destination state, when recorded). No-arrival
#'     maps have zero counts and NA ages.
#'   \item `counts`: age-bin by simulation counts, including zero-event maps.
#'   \item `summary`: bin counts (mean, median, lower, upper), `p.any` (fraction
#'     of maps with at least one event), and `prob.lower`, `prob.upper`.
#'   \item `overall`: probability of any matching event anywhere on the tree,
#'     its interval, and numbers of maps with and without an event.
#'   \item `count.distribution`: empirical total-event count probabilities.
#'   \item `age.summary`: first/last-arrival age medians and quantile intervals,
#'     conditional on at least one event, plus pooled-event age quantiles.
#'     First means oldest in forward history. First/last give each qualifying
#'     map one observation; pooled ages weight maps by their event counts.
#'   \item `event.age.distribution`: age-bin probabilities for a randomly
#'     selected event pooled across all maps (not probabilities of any event).
#'   \item `first.age.distribution`: bin probabilities for first arrival,
#'     conditional on at least one event; each qualifying map contributes once.
#' }
#' @details
#' Anagenetic events are read from within-edge maps. Cladogenetic changes
#' compare the recorded parental state with the initial state of each retained
#' daughter. Removed daughters never contribute. Two qualifying daughters at
#' the same split count as two lineage transitions; a shared ancestral branch
#' is counted once, not once per descendant tip. Root speciation is included
#' when its retained daughter changes from the recorded root state. Merely
#' starting in a destination state does not imply a migration event. If the
#' root already occupies the destination, its initial colonization can predate
#' the tree: first/last ages here describe recorded arrivals within the tree,
#' not that unobserved initial colonization.
#'
#' Count and age intervals are pointwise empirical simulation quantiles
#' (type 1 by default), not confidence intervals for their means. Probability intervals
#' are exact binomial Monte Carlo confidence intervals across independent maps;
#' they measure finite-simulation precision, not biological or parameter
#' uncertainty. All outputs describe the distribution of the supplied maps;
#' ASR-weighted maps are not exact ClaSSE posterior draws.
#'
#' Arrival means a state transition, not a directly observed dispersal event.
#' To include all arrivals onto S, supply every S-containing state in `to`,
#' not only S, A.S and S.R if other such states exist. `from` and `to` define
#' a cross-product of allowed directed changes. Events outside `breaks` are
#' reported in `per.map$n.outside` and retained in overall/age summaries, but
#' not in bins (a warning is issued). If no map has an arrival, conditional
#' age summaries and first-arrival probabilities are NA, not zero ages.
#' @param quantile.type Quantile algorithm (integer 1 to 9) for count and age
#'   summaries. Default 1 selects observed values; 7 interpolates. Does not
#'   change binomial probability intervals.
#' @export
transitions.through.time.simmap <- function(sims, breaks, to, from = NULL,
    event.types = c("anagenetic", "cladogenetic"), level = .95,
    quantile.type = 1L) {
  .classe_td_check_quantile_type(quantile.type)
  if (!is.null(sims$tree)) sims <- list(sims)
  if (!is.list(sims) || !length(sims)) stop("sims must contain at least one map")
  if (!is.numeric(breaks) || length(breaks)<2L || any(!is.finite(breaks)) ||
      any(breaks<0) || is.unsorted(breaks, strictly=TRUE))
    stop("breaks must contain at least two increasing nonnegative ages")
  valid.states <- function(x) is.character(x) && length(x)>0L && !anyNA(x) && !anyDuplicated(x)
  if (!valid.states(to) || (!is.null(from) && !valid.states(from)))
    stop("from and to must contain unique state labels")
  if (length(level)!=1L || !is.finite(level) || level<=0 || level>=1)
    stop("level must be between zero and one")
  event.types <- match.arg(event.types, c("anagenetic", "cladogenetic"), several.ok=TRUE)
  nb <- length(breaks)-1L
  n <- length(sims)
  counts <- matrix(0L, nb, n, dimnames=list(NULL, paste0("sim",seq_len(n))))
  records <- vector("list",n)
  per.map <- data.frame(simulation=seq_len(n), n.transitions=0L,
    n.anagenetic=0L, n.cladogenetic=0L, n.outside=0L,
    first.age=NA_real_, last.age=NA_real_, root.in.destination=NA)
  bin <- function(ages) as.integer(cut(ages, breaks, right=FALSE, include.lowest=TRUE))
  for (i in seq_len(n)) {
    labels <- colnames(sims[[i]]$tree$mapped.edge)
    if (!all(c(to,from) %in% labels)) stop("Unknown transition state in simulation ",i)
    source <- if (is.null(from)) setdiff(labels,to) else from
    e <- .classe_td_transition_events(sims[[i]], event.types)
    root <- setdiff(sims[[i]]$lineages$parent,sims[[i]]$lineages$child)
    root.event <- sims[[i]]$events
    if (!is.null(root.event$from.label)) {
      r <- which(root.event$node %in% root & root.event$event=="observed_speciation")
      if (length(r) && !is.na(root.event$from.label[r[1]]))
        per.map$root.in.destination[i] <- root.event$from.label[r[1]] %in% to
    }
    e <- e[e$from %in% source & e$to %in% to, , drop=FALSE]
    e$simulation <- rep.int(i,nrow(e))
    records[[i]] <- e
    b <- bin(e$age)
    counts[,i] <- tabulate(b,nbins=nb)
    per.map$n.transitions[i] <- nrow(e)
    per.map$n.anagenetic[i] <- sum(e$event.type=="anagenetic")
    per.map$n.cladogenetic[i] <- sum(e$event.type=="cladogenetic")
    per.map$n.outside[i] <- sum(is.na(b))
    if (nrow(e)) {
      per.map$first.age[i] <- max(e$age)
      per.map$last.age[i] <- min(e$age)
    }
  }
  if (any(per.map$n.outside)) warning("Some transition ages fall outside breaks; see per.map$n.outside")
  probs <- c((1-level)/2, 1-(1-level)/2)
  probability <- t(vapply(rowSums(counts>0), .classe_td_probability, numeric(3), n=n, level=level))
  summary <- data.frame(younger=head(breaks,-1), older=tail(breaks,-1),
    mean=rowMeans(counts), median=apply(counts,1,stats::median),
    lower=apply(counts,1,stats::quantile,probs=probs[1],type=quantile.type),
    upper=apply(counts,1,stats::quantile,probs=probs[2],type=quantile.type),
    p.any=probability[,1],prob.lower=probability[,2],prob.upper=probability[,3], n.maps=n)
  has <- per.map$n.transitions>0
  overall <- as.data.frame(as.list(.classe_td_probability(sum(has),n,level)))
  overall$n.maps <- n
  overall$n.with.event <- sum(has)
  overall$n.without.event <- sum(!has)
  overall$n.root.in.destination <- sum(per.map$root.in.destination,na.rm=TRUE)
  events <- do.call(rbind,records)
  age.summary <- do.call(rbind,lapply(c("first.age","last.age","all.events"), function(key) {
    v <- if(key=="all.events") events$age else per.map[[key]][has]
    q <- if(length(v)) stats::quantile(v,c(probs[1],.5,probs[2]),type=quantile.type) else rep(NA_real_,3)
    data.frame(which=key,median=q[2],lower=q[1],upper=q[3],n.maps=sum(has),
               n.ages=length(v),row.names=NULL)
  }))
  first.n <- tabulate(bin(per.map$first.age[has]),nbins=nb)
  first.p <- t(vapply(first.n,.classe_td_probability,numeric(3),n=sum(has),level=level))
  first <- data.frame(younger=head(breaks,-1),older=tail(breaks,-1),
    probability=first.p[,1],lower=first.p[,2],upper=first.p[,3],n.maps=sum(has))
  total <- tabulate(per.map$n.transitions+1L,nbins=max(per.map$n.transitions)+1L)
  event.age <- data.frame(younger=head(breaks,-1),older=tail(breaks,-1),
    n.events=rowSums(counts),probability=if(nrow(events)) rowSums(counts)/nrow(events) else NA_real_)
  structure(list(events=events,per.map=per.map,counts=counts,
    summary=summary,overall=overall,age.summary=age.summary,
    first.age.distribution=first,event.age.distribution=event.age,
    count.distribution=data.frame(n.transitions=seq_along(total)-1L,n.maps=total,probability=total/n),
    breaks=breaks,to=to,from=from,event.types=event.types,level=level,
    quantile.type=quantile.type),
    class="classe_td_transitions")
}

#' Plot transition probabilities or counts through time
#'
#' @param x Result of [transitions.through.time.simmap()].
#' @param type `"probability"` shows probability of any arrival per bin;
#'   `"counts"` shows median counts with simulation intervals;
#'   `"first.age"` shows conditional first-arrival bin probabilities.
#' @param col Plot colour.
#' @param xlab,ylab Axis labels.
#' @param ... Further arguments to the initial base plot.
#' @return Invisibly returns `x`.
#' @details Probability error bars are Monte Carlo confidence intervals;
#'   count error bars are map-distribution quantiles. Bins are not densities.
#' @export
plot.classe_td_transitions <- function(x,type=c("probability","counts","first.age"),
    col="#0072B2",xlab="Age before present",ylab=NULL,...) {
  type <- match.arg(type)
  s <- if(type=="first.age") x$first.age.distribution else x$summary
  if(type=="probability") {
    y <- s$p.any; lo <- s$prob.lower; hi <- s$prob.upper
  } else if(type=="counts") {
    y <- s$median; lo <- s$lower; hi <- s$upper
  } else { y <- s$probability; lo <- s$lower; hi <- s$upper }
  if(is.null(ylab)) ylab <- switch(type, probability="P(at least one arrival in bin)",
    counts="Number of transitions",first.age="P(first arrival in bin | any arrival)")
  top <- if(type=="counts") max(1,hi,y,na.rm=TRUE) else 1
  graphics::plot(NA,xlim=rev(range(x$breaks)),ylim=c(0,top),xlab=xlab,ylab=ylab,...)
  ok <- is.finite(y)
  mid <- (s$younger+s$older)/2
  if (any(ok)) {
    graphics::rect(s$younger[ok],0,s$older[ok],y[ok],
      col=grDevices::adjustcolor(col,.2),border=NA)
    graphics::segments(mid[ok],lo[ok],mid[ok],hi[ok],col=col)
    graphics::points(mid[ok],y[ok],pch=16,col=col)
  } else graphics::text(mean(range(x$breaks)),.5,"No arrivals: age distribution undefined")
  invisible(x)
}
