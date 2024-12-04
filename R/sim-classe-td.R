
table2tree <- function(info) {

  # Epoch 1
  tb.e1 <-  attr(info, "info.epoch1")
  phy.e1 <- diversitree:::me.to.ape.bisse(tb.e1[-1,], tb.e1$state[1])
  phy.e1 <- diversitree::prune(phy.e1)

  # all tree
  phy <- diversitree:::me.to.ape.bisse(info[-1,], info$state[1])
  phy <- prune(phy)
  # add info
  phy$epoch1.extant.tree.depth <- max(ape::node.depth.edgelength(phy.e1))[1]
  phy$epoch12.tree.depth <- max(ape::node.depth.edgelength(phy))[1]
  phy$epoch1.ntip <- ape::Ntip(phy.e1)
  phy$t.epoch1 <- attr(info, "t.epoch1")
  phy$t.total <- attr(info, "t.total")

  return(phy)
}

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

  #-------
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
  # if there no extant sp in epoch1 return Null
  tb.e1 <- diversitree:::me.to.ape.bisse(info1[-1,], info1$state[1])
  if (is.null(tb.e1))
    return(NULL)


  #--------------------
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

  # hist <- as.data.frame(do.call(rbind, hist))
  # if ( nrow(hist) == 0 )
  #   hist <- as.data.frame(matrix(NA, 0, 4))
  # names(hist) <- c("idx", "t", "from", "to")
  # hist$x0 <- info$start[match(hist$idx, info$idx)]
  # hist$tc <- hist$t - hist$x0

  #attr(info, "t") <- t
  attr(info, "info.epoch1") <- info1
  attr(info, "t.epoch1") <- t.epoch1
  attr(info, "t.total") <- t
  #attr(info, "hist") <- hist
  info
}





