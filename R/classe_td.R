
# Episodic ClaSSE.td model created using musse.td from diversitree package



# adopted from diversitree:::make.pars.musse.td
make.pars.classe.td <- function (n.epoch, k)
{
  # number of par in one epoch
  # np.in <- k * (k + 1)
  np.in <- (k^2+k^3)/2 + k + (k^2-k)

  # number of par in one epoch + time row (1 element) + N of negative elements in Q
  # np.out <- k * (k + 2) + 1
  np.out <- np.in + 1 + k

  # total number of pars: per_epoch * n of epach + 1
  npar <- (n.epoch - 1) + (np.in * n.epoch)

  i.t <- seq_len(n.epoch - 1)
  i.p <- n.epoch:npar

  #f.pars <- diversitree:::make.pars.musse(k)
  f.pars <- diversitree:::make.pars.classe(k)
  function(pars) {
    if (length(pars) != npar)
      stop(sprintf("Invalid length parameters (expected %d)",
                   npar))
    pars2 <- matrix(NA, n.epoch, np.out)
    pars2[, 1] <- c(pars[i.t], Inf)
    tmp <- matrix(pars[i.p], n.epoch, np.in, TRUE)
    for (i in seq_len(n.epoch)) pars2[i, -1] <- f.pars(tmp[i,
    ])
    pars2
  }
}


# diversitree:::make.classe.info
# function (k, phy)
# {
#   list(name = "classe", name.pretty = "ClaSSE", np = as.integer((k + 3) * k * k/2 + k),
#        argnames = default.argnames.classe(k),
#        ny = as.integer(2 * k), k = as.integer(k), idx.e = as.integer(1:k),
#        idx.d = as.integer((k + 1):(2 * k)), derivs = derivs.classe,
#        phy = phy, ml.default = "subplex", mcmc.lowerzero = TRUE,
#        doc = NULL, reference = c("Goldberg (submitted)"))
# }

# diversitree:::update.info.td
# info <- cache$info

#' @method update info.classe_td
#' @export
update.info.classe_td <- function (info, k, n.epoch)
{
  n.epoch <- diversitree:::check.n.epoch(n.epoch)

  argnames_classe = diversitree:::default.argnames.classe(k) #length(argnames_classe)
  argnames.td <- c(sprintf("t.%d", seq_len(n.epoch - 1)), diversitree:::argnames_twopart(argnames_classe, n.epoch))

  info$name = "classe"

  #info$np = as.integer((k + 3) * k * k/2 + k)
  # number of par in one epoch +  N of negative elements in Q
  info$np = (k^2+k^3)/2 + 2*k + (k^2-k)

  info$time.chunks <- TRUE
  info$argnames <- argnames.td

  info$derivs = diversitree:::derivs.classe
  info$ml.default = "subplex"
  info$mcmc.lowerzero = TRUE
  info$reference = c("Sergei Tarasov")

  #info$name.ode <- info$name
  info$name.pretty <- sprintf("%s (time-chunks)", 'ClaSSE')
  info$n.epoch <- n.epoch
  info
}

# diversitree:::make.cache.classe
# function (tree, states, k, sampling.f = NULL, strict = TRUE)
# {
#   if (k > 31)
#     stop("No more than 31 states allowed.  Increase in classe-eqs.c.")
#   cache <- diversitree:::make.cache.musse(tree, states, k, sampling.f, strict)
#   cache$info <- diversitree:::make.info.classe(k, tree)
#   cache
# }

make.cache.classe.td <- function (tree, states, k, n.epoch, sampling.f, strict)
{
  if (k > 31)
    stop("No more than 31 states allowed.  Increase in classe-eqs.c.")
  cache <- diversitree:::make.cache.musse(tree, states, k, sampling.f, strict)

  # cache$info <- diversitree:::update.info.td(cache$info, n.epoch)
  cache$info <- update.info.classe_td(cache$info, k, n.epoch)
  cache
}

# make.cache.classe.td <- function (tree, states, k, n.epoch, sampling.f, strict)
# {
#   if (k > 31)
#     stop("No more than 31 states allowed.  Increase in classe-eqs.c.")
#
#   #--- from make.cache.musse()
#   #cache <- diversitree:::make.cache.musse(tree, states, k, sampling.f, strict)
#   if (strict)
#     tree <- diversitree:::check.tree(tree)
#   states <- diversitree:::check.states(tree, states, strict=strict, strict.vals=1:k)
#   cache <- diversitree:::make.cache(tree)
#   cache$info <- diversitree:::make.info.musse(k, tree)
#   cache$states <- states
#   cache$sampling.f <- diversitree:::check.sampling.f(sampling.f, k)
#   cache$y <- diversitree:::initial.tip.xxsse(cache)
#   #---
#
#   # cache$info <- diversitree:::update.info.td(cache$info, n.epoch)
#   cache$info <- update.info.classe_td(cache$info, k, n.epoch)
#   cache
# }

# tree        = readRDS('my-test/data/phy.rds')
# states      = phy$tip.state+1
# k           = 3
# n.epoch     = 2
# sampling.f  = NULL
# strict      = TRUE
# control     = list()

#' Episodic ClaSSE model where different time epochs may have different parameters
#'
#' @param tree An ultrametric bifurcating phylogenetic tree, in ape “phylo” format.
#' @param states A vector of character states, from 1 to Inf. This vector must have names that correspond to the tip labels in the phylogenetic tree (tree$tip.label).
#' @param k The number of states. The maximum now is 31, but that can easily be increased if necessary.
#' @param n.epoch Number of epochs. 1 corresponds to plain ClaSSE, so this will generally be an integer at least 2.
#' @param sampling.f sampling fraction
#' @param strict The states vector is always checked to make sure that the values are integers on 1:k.
#' @param control List of control parameters for the ODE solver. See details in diversitree:::make.bisse.
#'
#' @return a function
#' @description
#' Prepare to run episodic ClaSSE on a phylogenetic tree and character states.
#' This function creates a likelihood function that can be used in ML or Bayesian inference.
#'
#' @details
#' This functions is based on diversitree::make.musse.td() and diversitree::make.classe().
#'
#'
#' @export
#'
#' @examples
#'
#' phy <- ape::read.tree(text = '(sp1:1.014858647,(((sp4:0.2614280769,sp5:0.2614280769)nd4:0.007869378189,sp3:0.2692974551)nd3:0.01878075832,sp2:0.2880782134)nd2:0.7267804339)nd1;')
#' states <- species_counts <- c(sp1 = 2, sp2 = 2, sp3 = 1, sp4 = 0, sp5 = 1)
#'
#' # classe model
#' lik.c <- diversitree::make.classe(phy, states+1, 3)
#' pars.c <- diversitree::starting.point.classe(phy, 3)
#' lik.c(pars.c)   # -8.367568
#'
#' # classe.td where two time bins have the same pars (=classe)
#' lik.ctd <- make.classe.td(phy, states + 1, k=3, n.epoch=2, control=list(backend = "gslode"))
#' par.ctd <- c(0.2, pars.c, pars.c)
#' names(par.ctd) <- diversitree::argnames(lik.ctd)
#' par.ctd
#' lik.ctd(par.ctd) # -8.367568
#'
#' # classe.td where two time bins have different pars
#' lik.ctd2 <- make.classe.td(phy, states + 1, k=3, n.epoch=2, control=list(backend = "gslode"))
#' par.ctd2 <- c(0.2, pars.c, pars.c)
#' names(par.ctd2) <- diversitree::argnames(lik.ctd)
#' par.ctd2['lambda222.2'] <- 0.4
#' par.ctd2['lambda222.1'] <- 0.1
#' lik.ctd(par.ctd2) # -8.268788
#'
make.classe.td <- function (tree, states, k, n.epoch, sampling.f = NULL, strict = TRUE, control = list())
{
  # cache <- diversitree:::make.cache.classe(tree, states, k, sampling.f, strict)
  # cache <- diversitree:::make.cache.musse.td(tree, states, k, n.epoch, sampling.f,strict)
  cache <- make.cache.classe.td(tree, states, k, n.epoch, sampling.f, strict)

  # initial.conditions
  initial.conditions <- diversitree:::make.initial.conditions.classe(k)
  all_branches <- diversitree:::make.all_branches.td.dtlik(cache, control, initial.conditions)

  #rootfunc <- make.rootfunc.td(cache, rootfunc.musse)
  rootfunc <- diversitree:::make.rootfunc.td(cache, diversitree:::rootfunc.classe)

  # f.pars <- diversitree:::make.pars.musse.td(n.epoch, k)
  # f.pars(pars.t)
  f.pars <- make.pars.classe.td(n.epoch, k)

  ll <- function(pars, condition.surv = TRUE, root = diversitree:::ROOT.OBS, root.p = NULL, intermediates = FALSE) {
    pars2 <- f.pars(pars)
    ans <- all_branches(pars2, intermediates)
    rootfunc(ans, pars2, condition.surv, root, root.p, intermediates)
  }

  class(ll) <- c("classe.td", "classe", "dtlik", "function")
  ll
}

