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

# Read in tree and Plot
library(divextra)
library(plyr)

phy <- readRDS('asr/tree.rds')
cols <- mapvalues(phy$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
plot(phy, label.offset = 2, cex=1, no.margin = F)
tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
nodelabels()
tiplabels()
edgelabels()
axisPhylo()

lik.td <-make.classe.td(phy, phy$tip.state_3, k=3, n.epoch=2, control=list(backend = "gslode"), strict=T)
par.td <- c(4, c(1:54)/100)
names(par.td) <- diversitree::argnames(lik.td)
lik.est <- lik.td(par.td, root=ROOT.FLAT, intermediates=T)
print(lik.est)

#----------
e <- environment(lik.td)
all_branches <- e$all_branches
rootfunc <- e$rootfunc

eb <- environment(e$all_branches)
branches <- eb$branches.td
initial.conditions <- eb$initial.conditions.td

f.pars <- e$f.pars
pars <- f.pars(par.td)

cache <- eb$cache
states.idx <- cache$info$idx.d

vals <- attr(lik.est, "intermediates")$vals
LQ <- attr(lik.est, "intermediates")$lq
rinit <- attr(lik.est, "intermediates")$init
base <- attr(lik.est, "intermediates")$base
res <- list(vals=vals, lq=LQ, init=rinit, base=base)

nodes = NULL
condition.surv = TRUE
root = ROOT.OBS
root.p = NULL

# diversitree:::do.asr.marginal.R

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

  # nd <- 9
  f <- function(nd) {
    ## Include current node but omit root:
    anc.nd <- c(nd, anc[[nd]])
    anc.nd <- anc.nd[-length(anc.nd)]
    p <- rep(NA, length(states.idx))

    # st=1
    for ( st in seq_along(states.idx) ) {
      #------
      relative.t.branch <- 1 # cannot be bigger branch length
      #abs.t.branch <- 1
      #------
      lq <- res$lq
      branch.init <- res$init
      branch.base <- res$base
      # branch.init[states.idx[-st],nd] <- 0
      y.in <- branch.init[,nd]

      # i=9
      for ( i in anc.nd ) {
        # there two ans (from 0 to relative.t.branch) and (from  relative.t.branch to len[i])
        # ans <- branches(y.in, len[i], pars, depth[i], i)
        ans_1 <- branches(y.in, relative.t.branch, pars, depth[i], i)
        ans_1.lq <- ans_1[[1]]
        ans_1.base <- ans_1[[2]]
        # y.in_2 is where Ancestral Probs are calculated
        y.in_2 <- c(ans_1.base)
        y.in_2[states.idx[-st]] <- 0
        ans_2 <- branches(y.in_2, len[i] - relative.t.branch, pars, depth[i] + relative.t.branch, i)
        ans_2.lq <- ans_2[[1]]
        ans_2.base <- ans_2[[2]]
        ans_total.lq <- ans_1.lq + ans_2.lq

        # lq[i] <- ans[[1]]
        # branch.base[,i] <- ans[[2]]
        lq[i] <- ans_total.lq
        branch.base[,i] <- ans_2.base
        #-----
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

