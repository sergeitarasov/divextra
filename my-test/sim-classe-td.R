library(yaml)
library(diversitree)
library(divextra)
source('R/utils-yaml.R')
source('R/utils-sse.R')


## GeoSSE equivalence
# Same tree simulated in ?make.geosse
pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
names(pars) <- diversitree:::default.argnames.geosse()
set.seed(5)
phy <- tree.geosse(pars, max.taxa=5, x0=0)
# saveRDS(phy, 'my-test/data/phy.rds')
# write.tree(phy, 'my-test/data/phy.tree')

root_state <- function() sample(0:2, 1) # root_state()
max.t1 = 4
max.t2 = 2
n.epoch <- 2
args.names <- diversitree:::default.argnames.geosse()
n.args <- n.epoch*length(args.names)

#--------
pars <- matrix(
  c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5,
    1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5
    ),
  n.epoch, n.args/n.epoch, byrow = T)
colnames(pars) <- args.names

# first epoch from the root
phy1 <- tree.geosse(pars[1,], max.t=max.t1, include.extinct=FALSE, x0= root_state() )
phy1$tip.state
plot(phy1)
nodelabels()
edgelabels()
phy1$edge.length

epoch1.ntips <- length(phy1$tip.state)
phy1.depth <- node.depth.edgelength(phy1)[1]





#-------------- Classe Function



## Parameters come in the order:
##   l_111, l112, ..., l_kkk, m_1, ..., m_k, qmat
## where qmat is the standard q matrix without the diagonals (same
## order as for make.classe).  As a result the parameters
## vector is (k*k*(k+1)/2) + (k) + (k*k-k) = (k+3)*k^2/2 elements long.

# pars.cl <- c(1:7)/10
# names(pars.cl) <- diversitree:::default.argnames.geosse()
# pars <- diversitree:::pars.ge.to.cl(pars.cl)

pars.cl <- c(1:27)/10
names(pars.cl) <- diversitree:::default.argnames.classe(3)
pars <- pars.cl

k=3
max.taxa= 5
max.t= Inf
#include.extinct=FALSE
x0=1
single.lineage=TRUE

tb <- diversitree:::make.tree.classe(pars, k=3, max.taxa=5, x0=1, single.lineage=TRUE)
tb <- diversitree:::make.tree.classe(pars, k=3, max.t=.5, x0=2, single.lineage=TRUE)
tb
str(tb)
attr(tb, "t")
attr(tb, "hist")

phy <- diversitree:::me.to.ape.bisse(tb[-1,], tb$state[1])
phy <- prune(phy)
plot(phy)
node.depth.edgelength(phy)[1]
node.depth.edgelength(phy)[1] + tb$len[1]
edgelabels(phy$edge.length)
phy$tip.state


#--- tb
pars.cl <- c(1:27)/10
names(pars.cl) <- diversitree:::default.argnames.classe(3)
pars <- pars.cl



#-- Args

pars.tb <- matrix(
  c(1:27/10,
    1:27/100
  ),
  2, 27, byrow = T)
colnames(pars.tb) <- diversitree:::default.argnames.classe(3)
pars.tb
max.t1=0.5
max.t2=1

tb <- make.tree.classe.tb(pars.tb, k=3, max.t1=1, max.t2=3, x0=1, single.lineage=TRUE)
tb
phy <- diversitree:::me.to.ape.bisse(tb[-1,], tb$state[1])
phy <- prune(phy)
plot.phylo(phy, show.tip.label =T)
phy$tip.state

max(node.depth.edgelength(phy))[1]
max(node.depth.edgelength(phy))[1]+ tb$len[1]
attr(tb, "t.epoch1")
attr(tb, "t.total")
attr(tb, "info.epoch1")

tb.e1 <-  attr(tb, "info.epoch1")
phy.e1 <- diversitree:::me.to.ape.bisse(tb.e1[-1,], tb.e1$state[1])
phy.e1 <- prune(phy.e1)
plot.phylo(phy.e1, show.tip.label =T)
max(node.depth.edgelength(phy.e1))[1]
max(node.depth.edgelength(phy.e1))[1]+ tb$len[1]
#--
tb <- make.tree.classe.tb(pars.tb, k=3, max.t1=1, max.t2=3, x0=1, single.lineage=TRUE)
tb
phy <- table2tree(tb)
plot(phy)

# tb.e1 <-  attr(tb, "info.epoch1")
# phy.e1 <- diversitree:::me.to.ape.bisse(tb.e1[-1,], tb.e1$state[1])
# phy.e1 <- diversitree::prune(phy.e1)


phy$t.epoch1
phy$t.total
phy$epoch1.extant.tree.depth
phy$epoch12.tree.depth
phy$epoch1.ntip
Ntip(phy)

table2tree <- function(info, include.extinct=FALSE) {

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

make.tree.classe.tb <- function(pars.tb, k, max.taxa=Inf, max.t1=Inf, max.t2=Inf, x0, single.lineage=TRUE) {
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

#-----------
pars.cl <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
names(pars.cl) <- diversitree:::default.argnames.geosse()
pars <- diversitree:::pars.ge.to.cl(pars.cl)

max.taxa= 5
max.t= Inf
include.extinct=FALSE
x0=1

phy.cs <-  tree.classe(pars, max.taxa=5, x0=1)
plot(phy.cs)

tree.classe <- function(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE, x0=NA) {
  L <- length(pars)
  k <- -1 + (L-1+sqrt(L*(L-2)))^(1/3) + (L-1-sqrt(L*(L-2)))^(1/3)
  # the solution above is exact, but it comes out "numeric"
  k <- as.integer(round(k))
  diversitree:::check.pars.classe(pars, k)

  if ( is.na(x0) )
    x0 <- sample(k, 1, FALSE, stationary.freq.classe(pars, k))
  if ( length(x0) != 1 || is.na(x0) || x0 < 1 || x0 > k )
    stop(paste("Invalid root state", x0))

  info <- diversitree:::make.tree.classe(pars, k, max.taxa, max.t, x0)
  #info

  phy <- diversitree:::me.to.ape.bisse(info[-1,], info$state[1])
  # plot(phy)
  # node.depth.edgelength(phy)[1]
  # edgelabels(phy$edge.length)
  # phy$tip.state

  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}



