# st <- asr.marginal(lik, pars)
#
# asr.marginal <- function(lik, pars, nodes=NULL, ...)
#   make.asr.marginal(lik)(pars, nodes, ...)
#
# make.asr.marginal <- function(lik, ...) {
#   UseMethod("make.asr.marginal")
# }
#
# diversitree:::make.asr.marginal.bisse(lik)(pars)


make.asr.marginal.bisse <- function(lik, ...) {

    e <- environment(lik)
    unresolved <- e$unresolved

    #-------
    do.asr <- diversitree:::make.do.asr.marginal(e$all_branches, e$rootfunc)

    all_branches <- e$all_branches
    rootfunc <- e$rootfunc
    #
    make.do.asr.marginal <- function (all_branches, rootfunc)
    {
      eb <- environment(all_branches)
      cache <- eb$cache
      states.idx <- cache$info$idx.d
      if (isTRUE(cache$info$partitioned)) {
        branches <- eb$branches.split
        initial.conditions <- eb$initial.conditions.split
      }
      else {
        branches <- eb$branches
        initial.conditions <- eb$initial.conditions
      }
      # preset <- NULL
      function(pars, nodes, preset, ...) {
        root.f <- function(res, pars) rootfunc(res, pars, ...)
        res <- all_branches(pars, TRUE, preset)
        do.asr.marginal.R(pars, cache, res, nodes, states.idx,
                          initial.conditions, branches, root.f)
      }
    }

    #---- asr
    asr <- function(pars, nodes = NULL, condition.surv = TRUE,
                    root = ROOT.FLAT, root.p = NULL, ...) {
      check.pars.bisse(pars)
      preset <- diversitree:::branches.unresolved.bisse(pars, unresolved)
      do.asr(pars, nodes, preset, condition.surv, root, root.p,
             intermediates = FALSE)
    }
    asr
}

nodes = NULL
condition.surv = TRUE
root = ROOT.OBS
root.p = NULL

diversitree:::do.asr.marginal.R
## Next, the utility functions for the different types of models This
## is to asr.marginal what all.branches is for the core models.  Here,
## the argument 'res' is the result of running all.branches
do.asr.marginal.R <- function(pars, cache, res, nodes, states.idx,
                              initial.conditions, branches, root,
                              ...) {
  ## Store these for easier calculation.
  children <- cache$children
  parent <- cache$parent
  len <- cache$len
  depth <- cache$depth
  root.idx <- cache$root
  anc <- cache$ancestors

  if ( is.null(nodes) )
    nodes <- root.idx:max(cache$order)
  else
    nodes <- nodes + cache$n.tip

  # nd <- 7
  f <- function(nd) {
    ## Include current node but omit root:
    anc.nd <- c(nd, anc[[nd]])
    anc.nd <- anc.nd[-length(anc.nd)]
    p <- rep(NA, length(states.idx))

    # st=1
    for ( st in seq_along(states.idx) ) {
      lq <- res$lq
      branch.init <- res$init
      branch.base <- res$base
      branch.init[states.idx[-st],nd] <- 0
      y.in <- branch.init[,nd]

      # i=7
      for ( i in anc.nd ) {
        ans <- branches(y.in, len[i], pars, depth[i], i)
        lq[i] <- ans[[1]]
        branch.base[,i] <- ans[[2]]
        j <- parent[i]
        y.in <- initial.conditions(branch.base[,children[j,]], pars,
                                   depth[j], j)
        branch.init[,j] <- y.in
      }

      # root <- root.f
      #root.f(list(vals=branch.init[,root.idx], lq=lq), pars)
      # !!!! root is mine here
      root <- function(res, pars) rootfunc(res, pars, root.p=NULL, root = ROOT.OBS, condition.surv=T, intermediates=F)
      #
      ans <- root(list(vals=branch.init[,root.idx], lq=lq), pars)

      ## explots IEEE arithmetic's exp(-Inf) == 0
      p[st] <- if ( is.na(ans) ) -Inf else ans
    }

    # x <- c(0.1, 0.2)
    # x/sum(x)
    # #
    # xl <- log(x)
    # xx <- exp(xl-max(xl))
    # xx/sum(xx)

    pp <- exp(p - max(p))
    pp/sum(pp)
  }

  matrix(unlist(lapply(nodes, f)), ncol=length(nodes))
}

