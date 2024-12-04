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



# second epoch from the root
phy.bind <- phy1
bind.tips.state <- c()
# tip=2
for (tip in 1:epoch1.ntips){
  print(tip)
  root.state <- phy1$tip.state[tip]
  phy2i <- tree.geosse(pars[2,], max.t=max.t2, include.extinct=FALSE, x0 = root.state )
  #plot(phy2i)
  if (!is.null(phy2i)){
    new.tips <- paste0(phy2i$tip.label, '.', tip)
    names(phy2i$tip.state) <- new.tips
    bind.tips.state <- c(bind.tips.state, phy2i$tip.state)
    phy2i$tip.label <- new.tips

    phy.bind <- bind.tree(phy.bind, phy2i, where = tip)
    #plot(phy.bind)
  }

}
plot(phy.bind)

str(phy2i)
str(phy.bind)
phy.bind$tip.state
phy.bind$tip.label
phy1 <- tree.geosse(pars[2,], max.t=max.t2, include.extinct=FALSE, x0= root_state() )
plot(phy1)


pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(3)
phy <- tree.bisse(pars, max.taxa=10, x0=0)
plot(phy)
node.depth.edgelength(phy)[1]
phy$tip.state




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

make.tree.classe <- function(pars, k, max.taxa=Inf, max.t=Inf, x0, single.lineage=TRUE) {
  # The other models don't require k, but this function is hidden away,
  # so no worry about passing k in rather than recomputing it.

  if ( x0 < 1 || x0 > k )
    stop("x0 must be an integer in [1,k]")

  # arrange the parameters in a list with elements:
  #   lambda = lambda_ijk array, mu = mu vector, q = q_ij array,
  #   nstates = number of states
  pars.list <- diversitree:::inflate.pars.classe(pars, k)

  # for drawing samples below, it's nicer to have 0 than NA for the
  # non-applicable speciation rates
  pars.list$lambda[which(is.na(pars.list$lambda))] <- 0

  # row i is all states != i, i.e., states that can be transitioned to
  to <- matrix(unlist(lapply(1:k, function(i) (1:k)[-i])), k, k-1, TRUE)

  # pars is a "k x k+1" matrix giving, for a lineage in state row i,
  # the rate at which speciation, extinction, or anagenetic transition
  # to each other state happens.  This approach loses speciation info
  # (retained in pars.list$lambda) and requires an extra sample() call
  # within the speciation "if" below, but it makes the indices less
  # heinous.

  # cols are: speciation, extincton, rest are changes to states
  pars <- cbind(rowSums(pars.list$lambda), pars.list$mu,
                matrix(pars[-seq_len(k*k*(k+1)/2+k)], k, k-1, TRUE))
  # r.i = total rate at which something happens to a lineage in state i
  r.i <- rowSums(pars)

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

  while ( n.taxa <= max.taxa && n.taxa > 0 ) {

    # When does an event happen?
    r.n <- r.i * n.i
    r.tot <- sum(r.n)
    dt <- rexp(1, r.tot)
    t <- t + dt

    # Stop if it happens too late.
    if ( t > max.t ) {
      dt <- dt - (t - max.t)
      len[lineages] <- len[lineages] + dt
      t <- max.t
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
    type <- sample(k+1, 1, FALSE, pars[state,])

    if ( type == 1 ) {                      # Speciation
      if ( n.taxa == max.taxa )
        break
      new.i <- length(extinct) + 1:2
      split[lineage] <- TRUE
      extinct[new.i] <- split[new.i] <- FALSE

      # get daughter states from indices of lamda_ijk
      lam <- pars.list$lambda[state,,]
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
      states[lineage] <- state.new <- to[state, type - 2]
      n.i[c(state.new, state)] <- n.i[c(state.new, state)] + c(1,-1)
      hist[[length(hist)+1]] <- c(lineage, t, state, state.new)
    }
  }

  info <- data.frame(idx=seq_along(extinct), len=len, parent=parent,
                     start=start, state=states, extinct=extinct,
                     split=split)

  hist <- as.data.frame(do.call(rbind, hist))
  if ( nrow(hist) == 0 )
    hist <- as.data.frame(matrix(NA, 0, 4))
  names(hist) <- c("idx", "t", "from", "to")
  hist$x0 <- info$start[match(hist$idx, info$idx)]
  hist$tc <- hist$t - hist$x0

  attr(info, "t") <- t
  attr(info, "hist") <- hist
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
  info

  phy <- diversitree:::me.to.ape.bisse(info[-1,], info$state[1])
  plot(phy)
  node.depth.edgelength(phy)[1]
  edgelabels(phy$edge.length)
  phy$tip.state

  if ( include.extinct || is.null(phy) )
    phy
  else
    prune(phy)
}



