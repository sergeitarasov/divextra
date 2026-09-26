

#' Converts GeoSSE parameters to ClaSSE
#'
#' @param pars.ge GeoSSE parameters
#'
#' @return vector
#' @description
#' It is based on diversitree:::pars.ge.to.cl(). Only works for three states
#' and time-homogeneous models. The state order is important: A, B, AB.
#'
#' @export
#'
#' @examples
#' ## "sA"  "sB" "sAB" "xA"  "xB"  "dA"  "dB"
#' pars.ge <- matrix(
#' c(0.1, 0.1, 0.1,  0.001, 0.001,  0.1, 0.1,
#'   0.3, 0.3, 0.3,  0.001, 0.001,  0.1, 0.1
#' ),
#' 2, 7, byrow = TRUE)
#' colnames(pars.ge) <- diversitree:::default.argnames.geosse()
#' pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))
#' print(pars.cl)
pars.geosse2classe <- function(pars.ge)
{
  if (is.null(names(pars.ge)))
    names(pars.ge) <- diversitree:::default.argnames.geosse()
  pars.cl <- rep(0, 27)
  names(pars.cl) <- diversitree:::default.argnames.classe(3)
  pars.cl['lambda111'] <- pars.cl['lambda313'] <- pars.ge['sA']
  pars.cl['lambda222'] <- pars.cl['lambda323'] <- pars.ge['sB']
  pars.cl['lambda312'] <-  pars.ge['sAB']
  pars.cl['mu1'] <- pars.cl['q32'] <- pars.ge['xA']
  pars.cl['mu2'] <- pars.cl['q31'] <- pars.ge['xB']
  pars.cl['q13'] <- pars.ge['dA']
  pars.cl['q23'] <- pars.ge['dB']
  pars.cl
}


.classe_td_history_table <- function(hist, info) {
  if (!length(hist)) {
    return(data.frame(
      idx = integer(), t = numeric(), from = integer(), to = integer(),
      x0 = numeric(), tc = numeric()
    ))
  }

  ans <- as.data.frame(do.call(rbind, hist))
  names(ans) <- c("idx", "t", "from", "to")
  ans$idx <- as.integer(ans$idx)
  ans$from <- as.integer(ans$from)
  ans$to <- as.integer(ans$to)
  ans$x0 <- info$start[match(ans$idx, info$idx)]
  ans$tc <- ans$t - ans$x0
  ans
}

.classe_td_join_maps <- function(first, second) {
  if (!length(first))
    return(second)
  if (!length(second))
    return(first)
  if (identical(names(first)[length(first)], names(second)[1])) {
    first[length(first)] <- first[length(first)] + second[1]
    if (length(second) > 1L)
      first <- c(first, second[-1L])
    first
  } else {
    c(first, second)
  }
}

.classe_td_drop_tip_simmap <- function(phy, tip) {
  old.n.tip <- ape::Ntip(phy)
  drop <- if (is.character(tip)) match(tip, phy$tip.label) else as.integer(tip)
  drop <- drop[!is.na(drop)]
  if (!length(drop))
    return(phy)
  if (old.n.tip - length(unique(drop)) < 2L)
    return(NULL)

  remove <- match(drop, phy$edge[, 2])
  keep <- setdiff(seq_len(nrow(phy$edge)), remove)
  phy$edge <- phy$edge[keep, , drop = FALSE]
  phy$edge.length <- phy$edge.length[keep]
  phy$maps <- phy$maps[keep]

  # Remove ancestral edges that no longer lead to a retained terminal.
  repeat {
    dead <- setdiff(phy$edge[, 2], phy$edge[, 1])
    dead <- dead[dead > old.n.tip]
    if (!length(dead))
      break
    remove <- match(dead, phy$edge[, 2])
    keep <- setdiff(seq_len(nrow(phy$edge)), remove)
    phy$edge <- phy$edge[keep, , drop = FALSE]
    phy$edge.length <- phy$edge.length[keep]
    phy$maps <- phy$maps[keep]
  }

  # Discard a leading stem, matching table2tree(), then concatenate every
  # remaining unary chain in root-to-tip map order.
  repeat {
    root <- setdiff(phy$edge[, 1], phy$edge[, 2])
    root.edges <- which(phy$edge[, 1] == root)
    if (length(root.edges) != 1L)
      break
    remove <- root.edges
    keep <- setdiff(seq_len(nrow(phy$edge)), remove)
    phy$edge <- phy$edge[keep, , drop = FALSE]
    phy$edge.length <- phy$edge.length[keep]
    phy$maps <- phy$maps[keep]
  }

  repeat {
    parents <- unique(phy$edge[, 1])
    out.degree <- vapply(parents, function(node) {
      sum(phy$edge[, 1] == node)
    }, integer(1))
    unary <- parents[out.degree == 1L & parents %in% phy$edge[, 2]]
    if (!length(unary))
      break
    node <- unary[1]
    incoming <- which(phy$edge[, 2] == node)
    outgoing <- which(phy$edge[, 1] == node)
    phy$edge[incoming, 2] <- phy$edge[outgoing, 2]
    phy$maps[[incoming]] <- .classe_td_join_maps(
      phy$maps[[incoming]], phy$maps[[outgoing]]
    )
    phy$edge.length[incoming] <- sum(phy$maps[[incoming]])
    keep <- setdiff(seq_len(nrow(phy$edge)), outgoing)
    phy$edge <- phy$edge[keep, , drop = FALSE]
    phy$edge.length <- phy$edge.length[keep]
    phy$maps <- phy$maps[keep]
  }

  terminal.old <- setdiff(phy$edge[, 2], phy$edge[, 1])
  phy$tip.label <- phy$tip.label[terminal.old]
  phy$edge[match(terminal.old, phy$edge[, 2]), 2] <- seq_along(terminal.old)
  internal.old <- sort(unique(phy$edge[, 1]))
  internal.new <- seq_along(internal.old) + length(terminal.old)
  for (z in seq_along(internal.old))
    phy$edge[phy$edge == internal.old[z]] <- internal.new[z]
  storage.mode(phy$edge) <- "integer"
  phy$Nnode <- length(internal.old)
  phy$node.label <- NULL
  phy$tip.state <- NULL
  phy$node.state <- NULL
  phy$edge.state <- NULL
  phy$mapped.edge <- NULL
  phy
}


#' Converts episodic ClaSSE.td simulations to a tree
#'
#' @param info a table from make.tree.classe.td()
#'
#' @return tree
#' @description
#' Converts an output (a table) from make.tree.classe.td() to a phylogenetic tree
#' that contains extant species only. It is based on diversitree:::me.to.ape.bisse().
#'
#' @export
#'
#' @examples
#'
#' ## "sA"  "sB" "sAB" "xA"  "xB"  "dA"  "dB"
#' pars.ge <- matrix(
#'   c(0.1, 0.1, 0.1,  0, 0,  0.1, 0.1,
#'     0.3, 0.3, 0.3,  0, 0,  0.1, 0.1
#'   ),
#'   2, 7, byrow = TRUE)
#' colnames(pars.ge) <- diversitree:::default.argnames.geosse()
#' pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))
#' print(pars.cl)
#'
#' set.seed(123)
#' tb <- make.tree.classe.td(pars.cl, k=3, max.t1=10, max.t2=15, x0=1, single.lineage=TRUE)
#' print(tb)
#' phy <- table2tree(tb)
#' plot(phy)
#'
table2tree <- function(info) {

  # Epoch 1
  tb.e1 <-  attr(info, "info.epoch1")
  phy.e1 <- diversitree:::me.to.ape.bisse(tb.e1[-1,], tb.e1$state[1])
  phy.e1 <- diversitree::prune(phy.e1)

  # all tree
  phy <- diversitree:::me.to.ape.bisse(info[-1,], info$state[1])
  phy <- diversitree::prune(phy)
  # add info
  phy$epoch1.extant.tree.depth <- max(ape::node.depth.edgelength(phy.e1))[1]
  phy$epoch12.tree.depth <- max(ape::node.depth.edgelength(phy))[1]
  phy$epoch1.ntip <- ape::Ntip(phy.e1)
  phy$t.epoch1 <- attr(info, "t.epoch1")
  phy$t.total <- attr(info, "t.total")
  #
  # In simulations epoch1 starts  at the root, while in inference it starts at the tips
  # Thus, for the inference the model switches regimes at
  phy$t.regime.change <- phy$epoch12.tree.depth - phy$epoch1.extant.tree.depth
  phy$sim.pars <- attr(info, "sim.pars") # parameters used in simulation

  return(phy)
}


#' Convert an episodic ClaSSE simulation table to a stochastic map
#'
#' @param info A table returned by [make.tree.classe.td()]. Simulations made
#'   with older versions of divextra do not contain the transition-history
#'   attribute and can only be converted to constant-state edge maps.
#' @param state.labels Optional character vector giving the plotted label for
#'   states `1, ..., k`. For the three-state GeoSSE parameterization, for
#'   example, use `c("A", "B", "AB")`.
#'
#' @return A tree of class `simmap` and `phylo`. Its `maps` element contains
#'   the ordered state durations along every retained edge and `mapped.edge`
#'   contains the total duration in each state.
#'
#' @description
#' Converts the complete event table from [make.tree.classe.td()] into a
#' `phytools`-compatible stochastic mapping tree. Extinct species are pruned in
#' the same way as in [table2tree()]. Anagenetic transitions on retained and
#' collapsed branches are preserved.
#'
#' @export
#'
#' @examples
#' pars.ge <- matrix(
#'   c(0.1, 0.1, 0.1, 0, 0, 0.1, 0.1,
#'     0.3, 0.3, 0.3, 0, 0, 0.1, 0.1),
#'   2, 7, byrow = TRUE
#' )
#' colnames(pars.ge) <- diversitree:::default.argnames.geosse()
#' pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))
#'
#' set.seed(123)
#' tb <- make.tree.classe.td(
#'   pars.cl, k = 3, max.t1 = 10, max.t2 = 15,
#'   x0 = 1, single.lineage = TRUE
#' )
#' simmap <- table2simmap(tb, state.labels = c("A", "B", "AB"))
#' if (requireNamespace("phytools", quietly = TRUE)) {
#'   colors <- c(A = "#D55E00", B = "#0072B2", AB = "#009E73")
#'   phytools::plotSimmap(simmap, colors = colors, ftype = "i")
#' }
table2simmap <- function(info, state.labels = NULL) {
  full <- diversitree:::me.to.ape.bisse(info[-1,], info$state[1])
  if (is.null(full))
    return(NULL)

  hist <- full$hist
  states <- c(
    unname(full$tip.state), unname(full$node.state),
    if (!is.null(hist) && nrow(hist)) c(hist$from, hist$to)
  )
  states <- as.integer(states[is.finite(states)])
  k <- attr(info, "k")
  if (is.null(k))
    k <- if (length(states)) max(states) else 0L
  k <- as.integer(k)
  if (k < 1L)
    stop("Could not determine the number of states")

  if (is.null(state.labels))
    state.labels <- as.character(seq_len(k))
  if (length(state.labels) != k || anyNA(state.labels) ||
      any(!nzchar(state.labels)) || anyDuplicated(state.labels)) {
    stop("state.labels must contain k unique, non-empty labels")
  }
  state.labels <- as.character(state.labels)

  if (is.null(hist)) {
    warning(
      "No transition history is stored in info; creating constant-state ",
      "edge maps. Regenerate the table with make.tree.classe.td() to retain ",
      "within-edge transitions."
    )
  }

  node.names <- c(full$tip.label, full$node.label)
  terminal.states <- c(unname(full$tip.state), unname(full$node.state))
  tol <- 1e-10

  maps <- lapply(seq_len(nrow(full$edge)), function(edge.index) {
    child <- full$edge[edge.index, 2]
    child.name <- node.names[child]
    edge.length <- full$edge.length[edge.index]
    edge.hist <- if (is.null(hist) || !nrow(hist)) {
      NULL
    } else {
      hist[!is.na(hist$name2) & hist$name2 == child.name, , drop = FALSE]
    }

    if (!is.null(edge.hist) && nrow(edge.hist)) {
      edge.hist <- edge.hist[order(edge.hist$tc), , drop = FALSE]
      if (any(edge.hist$tc < -tol | edge.hist$tc > edge.length + tol))
        stop("Stored transition time lies outside its retained edge")
      transition.times <- pmin(edge.length, pmax(0, edge.hist$tc))
      map.states <- c(edge.hist$from[1], edge.hist$to)
    } else {
      transition.times <- numeric()
      map.states <- terminal.states[child]
    }

    durations <- diff(c(0, transition.times, edge.length))
    positive <- durations > tol
    durations <- durations[positive]
    map.states <- map.states[positive]
    if (!length(durations)) {
      durations <- edge.length
      map.states <- terminal.states[child]
    }

    map <- numeric()
    for (z in seq_along(durations)) {
      label <- state.labels[map.states[z]]
      if (length(map) && identical(names(map)[length(map)], label)) {
        map[length(map)] <- map[length(map)] + durations[z]
      } else {
        map <- c(map, stats::setNames(durations[z], label))
      }
    }
    map
  })

  full$maps <- maps
  full$mapped.edge <- matrix(
    0, nrow = nrow(full$edge), ncol = k,
    dimnames = list(
      paste(full$edge[, 1], full$edge[, 2], sep = ","), state.labels
    )
  )
  for (edge.index in seq_along(full$maps)) {
    totals <- tapply(
      as.numeric(full$maps[[edge.index]]),
      names(full$maps[[edge.index]]), sum
    )
    full$mapped.edge[edge.index, names(totals)] <- totals
  }
  class(full) <- c("simmap", "phylo")

  extinct <- full$tip.label[startsWith(full$tip.label, "ex")]
  if (length(extinct)) {
    phy <- .classe_td_drop_tip_simmap(full, extinct)
    if (is.null(phy))
      return(NULL)
  } else {
    phy <- full
  }

  mapped.edge <- matrix(
    0, nrow = nrow(phy$edge), ncol = k,
    dimnames = list(
      paste(phy$edge[, 1], phy$edge[, 2], sep = ","), state.labels
    )
  )
  for (edge.index in seq_along(phy$maps)) {
    totals <- tapply(
      as.numeric(phy$maps[[edge.index]]), names(phy$maps[[edge.index]]), sum
    )
    mapped.edge[edge.index, names(totals)] <- totals
  }

  phy$mapped.edge <- mapped.edge
  phy$state.labels <- state.labels
  phy$tip.state <- full$tip.state[phy$tip.label]
  phy$edge.state <- match(
    vapply(phy$maps, function(map) names(map)[length(map)], character(1)),
    state.labels
  )

  # Pruning a simmap may collapse a chain of original lineages. The maps retain
  # those segments, including instantaneous cladogenetic state changes.
  phy$node.label <- sprintf("nd%d", seq_len(phy$Nnode))
  phy$node.state <- stats::setNames(rep(NA_integer_, phy$Nnode), phy$node.label)
  root <- setdiff(phy$edge[, 1], phy$edge[, 2])
  phy$node.state[root - ape::Ntip(phy)] <- info$state[1]
  internal.edges <- which(phy$edge[, 2] > ape::Ntip(phy))
  phy$node.state[phy$edge[internal.edges, 2] - ape::Ntip(phy)] <-
    phy$edge.state[internal.edges]
  epoch1 <- attr(info, "info.epoch1")
  phy$epoch1.extant.tree.depth <- attr(info, "t.epoch1") - epoch1$len[1]
  phy$epoch12.tree.depth <- attr(info, "t.total") - info$len[1]
  phy$epoch1.ntip <- sum(!epoch1$split & !epoch1$extinct)
  phy$t.epoch1 <- attr(info, "t.epoch1")
  phy$t.total <- attr(info, "t.total")
  phy$t.regime.change <- phy$t.total - phy$t.epoch1
  phy$sim.pars <- attr(info, "sim.pars")
  phy$simulation.table <- info
  phy$transition.history <- attr(info, "hist")
  class(phy) <- c("simmap", "phylo")
  phy
}





#' Simulate episodic ClaSSE tree and character states
#'
#' @param pars.tb parameters, a matrix with two rows with each row correponding to one epoch
#' @param k number of states
#' @param max.taxa Maximum number of taxa to include in the tree. If Inf, then the tree will be evolved until max.t time has passed.
#' @param max.t1 	Maximum tree length to evolve during epoch 1.
#' @param max.t2 Maximum tree length to evolve in total (i.e., epoch 1 + epoch 2)
#' @param x0 Initial character state at the root (1, .., k)
#' @param single.lineage Start simulation with a single lineage?
#'
#' @return data frame
#' @description
#' It is based on diversitree:::make.tree.classe() but simulates an episodic (time-dependent) ClaSSE tree
#' where the tree is split into two time epochs, each with its own parameters.
#'
#' @details
#' The time epochs start from the root, so the simulation begins with epoch 1 and ends in epoch 2.
#' The maximum number of epochs is two.' Note, in the ML calculation with make.classe.td(),
#' the epoch order is reversed -- epoch 1 starts at tips.
#'
#' This function simulates the tree using the code from diversitree:::make.tree.classe() for epoch one defined by
#' max.t1. Next, the tree from epoch one is used to simulate the final tree with the code from diversitree:::make.tree.classe() as well.
#'
#' The simulations are tested only with max.t1 and max.t2 arguments. Use max.taxa at your own risk.
#' The argument single.lineage should be always set to TRUE.
#'
#'
#' The returned data frame stores the first-epoch table, epoch times, simulation
#' parameters, number of states, and the complete anagenetic transition history
#' as attributes. The history columns are lineage index, absolute event time,
#' source state, destination state, lineage start time, and within-lineage time.
#'
#' @export
#'
#' @examples
#'
#' ## "sA"  "sB" "sAB" "xA"  "xB"  "dA"  "dB"
#' pars.ge <- matrix(
#'   c(0.1, 0.1, 0.1,  0, 0,  0.1, 0.1,
#'     0.3, 0.3, 0.3,  0, 0,  0.1, 0.1
#'   ),
#'   2, 7, byrow = TRUE)
#' colnames(pars.ge) <- diversitree:::default.argnames.geosse()
#' pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))
#' print(pars.cl)
#'
#' set.seed(123)
#' tb <- make.tree.classe.td(pars.cl, k=3, max.t1=10, max.t2=15, x0=1, single.lineage=TRUE)
#' print(tb)
#' phy <- table2tree(tb)
#' plot(phy)
#'
#' # the tree from epoch 1 only
#' tb.e1 <-  attr(tb, "info.epoch1")
#' phy.e1 <- diversitree:::me.to.ape.bisse(tb.e1[-1,], tb.e1$state[1])
#' phy.e1 <- diversitree::prune(phy.e1) # prune extinct
#' plot(phy.e1)
#'
#' # the final tree (after epoch 2)
#' phy$t.epoch1 # time of epoch 1 (=max.t1)
#' phy$t.total # total time (=max.t2)
#' phy$epoch1.extant.tree.depth # tree depth from epoch 1
#' phy$epoch12.tree.depth # total tree depth
#' phy$epoch1.ntip # N of tips in epoch 1 only
#' ape::Ntip(phy) # total N of tips
#'
#' # The model changes regimes at (if counting from the tips for the infference):
#' phy$t.regime.change # := phy$epoch12.tree.depth - phy$epoch1.extant.tree.depth
make.tree.classe.td <- function(pars.tb, k, max.taxa=Inf, max.t1=Inf, max.t2=Inf, x0, single.lineage=TRUE) {
  # The other models don't require k, but this function is hidden away,
  # so no worry about passing k in rather than recomputing it.

  pars1 <- pars.tb[1,]

  if ( x0 < 1 || x0 > k )
    stop("x0 must be an integer in [1,k]")

  # arrange the parameters in a list with elements:
  #   lambda = lambda_ijk array, mu = mu vector, q = q_ij array,
  #   nstates = number of states
  #pars.list <- diversitree:::inflate.pars.classe(pars, k)
  pars1.list <- diversitree:::inflate.pars.classe(pars1, k)

  # for drawing samples below, it's nicer to have 0 than NA for the
  # non-applicable speciation rates
  #pars.list$lambda[which(is.na(pars.list$lambda))] <- 0
  pars1.list$lambda[which(is.na(pars1.list$lambda))] <- 0

  # row i is all states != i, i.e., states that can be transitioned to
  #to <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])), k, k-1, TRUE)
  to1 <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])), k, k-1, TRUE)

  # pars is a "k x k+1" matrix giving, for a lineage in state row i,
  # the rate at which speciation, extinction, or anagenetic transition
  # to each other state happens.  This approach loses speciation info
  # (retained in pars.list$lambda) and requires an extra sample() call
  # within the speciation "if" below, but it makes the indices less
  # heinous.

  # cols are: speciation, extincton, rest are changes to states
  pars1 <- cbind(rowSums(pars1.list$lambda), pars1.list$mu,
                matrix(pars1[-seq_len(k*k*(k+1)/2+k)], k, k-1, TRUE))
  # r.i = total rate at which something happens to a lineage in state i
  r.i.1 <- rowSums(pars1)

  #------- Initialization
  extinct <- FALSE
  split   <- FALSE
  parent <- 0
  n.i <- rep(0, k)    # number of lineages in state i at this time
  len <- 0            # branch lengths
  t <- 0              # time elapsed
  hist <- list()      # history of transitions

  if ( single.lineage ) {
    states <- x0
    n.taxa <- lineages <- n.i[x0] <- 1
    start <- 0
  } else {
    ##states <- rep(x0, 2)
    ##n.taxa <- lineages <- n.i[x0] <- 2
    stop("Nope.")
  }

  # -------- Epoch 1
  while ( n.taxa <= max.taxa && n.taxa > 0 ) {

    # When does an event happen?
    r.n <- r.i.1 * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    # Stop if it happens too late.
    if ( t > max.t1 ) {
      dt <- dt - (t - max.t1)
      len[lineages] <- len[lineages] + dt
      t <- max.t1
      break
    }

    len[lineages] <- len[lineages] + dt

    # What state does the event happen to?
    state <- sample(k, 1, FALSE, r.n/r.tot)

    # What lineage with that state gets the event?
    j <- sample(n.i[state], 1)
    lineage <- lineages[states[lineages] == state][j]

    # What event happens?  1 = speciation, 2 = extinction,
    #    type>2 = transition (type & to provide new state)
    type <- sample(k+1, 1, FALSE, pars1[state,])

    if ( type == 1 ) {                      # Speciation
      if ( n.taxa == max.taxa )
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE

      # get daughter states from indices of lamda_ijk
      lam <- pars1.list$lambda[state,,]
      j <- sample(k*k, 1, FALSE, lam)
      s.daught <- c((j-1) %% k + 1, (j-1) %/% k + 1)
      states[new.i] <- s.daught
      n.i[state] <- n.i[state] - 1
      n.i[s.daught[1]] <- n.i[s.daught[1]] + 1
      n.i[s.daught[2]] <- n.i[s.daught[2]] + 1

      parent[new.i] <- lineage
      start[new.i] <- t
      len[new.i] <- 0
      n.taxa <- n.taxa + 1
      lineages <- which(!split & !extinct)

    } else if ( type == 2 ) {               # Extinction
      extinct[lineage] <- TRUE
      lineages <- which(!split & !extinct)
      n.i[state] <- n.i[state] - 1
      n.taxa <- n.taxa - 1

    } else {                                # Transition (anagenetic)
      states[lineage] <- state.new <- to1[state, type - 2]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1,-1)
      hist[[length(hist)+1]] <- c(lineage, t, state, state.new)
    }
  }

  info1 <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     start=start, state=states, extinct=extinct,
                     split=split)
  attr(info1, "hist") <- .classe_td_history_table(hist, info1)
  attr(info1, "k") <- k
  # if there no extant sp in epoch1 return Null
  tb.e1 <- diversitree:::me.to.ape.bisse(info1[-1,], info1$state[1])
  if (is.null(tb.e1))
    return(NULL)


  #-------------------- setup for Epoch 2
  t.epoch1 <-  t

  # setup for Epoch 2
  pars2 <- pars.tb[2,]
  pars2.list <- diversitree:::inflate.pars.classe(pars2, k)
  pars2.list$lambda[which(is.na(pars2.list$lambda))] <- 0
  to2 <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])), k, k-1, TRUE)
  pars2 <- cbind(rowSums(pars2.list$lambda), pars2.list$mu,
                 matrix(pars2[-seq_len(k*k*(k+1)/2+k)], k, k-1, TRUE))
  r.i.2 <- rowSums(pars2)

  # -------- Epoch 2
  while ( n.taxa <= max.taxa && n.taxa > 0 ) {

    # When does an event happen?
    r.n <- r.i.2 * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    # Stop if it happens too late.
    if ( t > max.t2 ) {
      dt <- dt - (t - max.t2)
      len[lineages] <- len[lineages] + dt
      t <- max.t2
      break
    }

    len[lineages] <- len[lineages] + dt

    # What state does the event happen to?
    state <- sample(k, 1, FALSE, r.n/r.tot)

    # What lineage with that state gets the event?
    j <- sample(n.i[state], 1)
    lineage <- lineages[states[lineages] == state][j]

    # What event happens?  1 = speciation, 2 = extinction,
    #    type>2 = transition (type & to provide new state)
    type <- sample(k+1, 1, FALSE, pars2[state,])

    if ( type == 1 ) {                      # Speciation
      if ( n.taxa == max.taxa )
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE

      # get daughter states from indices of lamda_ijk
      lam <- pars2.list$lambda[state,,]
      j <- sample(k*k, 1, FALSE, lam)
      s.daught <- c((j-1) %% k + 1, (j-1) %/% k + 1)
      states[new.i] <- s.daught
      n.i[state] <- n.i[state] - 1
      n.i[s.daught[1]] <- n.i[s.daught[1]] + 1
      n.i[s.daught[2]] <- n.i[s.daught[2]] + 1

      parent[new.i] <- lineage
      start[new.i] <- t
      len[new.i] <- 0
      n.taxa <- n.taxa + 1
      lineages <- which(!split & !extinct)

    } else if ( type == 2 ) {               # Extinction
      extinct[lineage] <- TRUE
      lineages <- which(!split & !extinct)
      n.i[state] <- n.i[state] - 1
      n.taxa <- n.taxa - 1

    } else {                                # Transition (anagenetic)
      states[lineage] <- state.new <- to2[state, type - 2]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1,-1)
      hist[[length(hist)+1]] <- c(lineage, t, state, state.new)
    }
  }

  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     start=start, state=states, extinct=extinct,
                     split=split)

  #attr(info, "t") <- t
  attr(info, "info.epoch1") <- info1
  attr(info, "t.epoch1") <- t.epoch1 # time of epoch 1 (=max.t1)
  attr(info, "t.total") <- t # total time (=max.t2)
  attr(info, "sim.pars") <- pars.tb # parameters used in the simulation
  attr(info, "hist") <- .classe_td_history_table(hist, info)
  attr(info, "k") <- k
  info
}
