

make.do.asr.marginal_divextra <- function(all_branches, rootfunc)
{
  eb <- environment(all_branches)
  cache <- eb$cache
  states.idx <- cache$info$idx.d
  if (isTRUE(cache$info$partitioned)) {
    branches <- eb$branches.split
    initial.conditions <- eb$initial.conditions.split
  }
  else {
    # branches <- eb$branches
    branches <- eb$branches.td
    #initial.conditions <- eb$initial.conditions
    initial.conditions <- eb$initial.conditions.td
  }
  function(pars, nodes, preset, ...) {
    root.f <- function(res, pars) rootfunc(res, pars, ...)
    res <- all_branches(pars, TRUE, preset)
    diversitree:::do.asr.marginal.R(pars, cache, res, nodes, states.idx,
                      initial.conditions, branches, root.f)
  }
}



make.asr.marginal.classe <- function (lik, ...)
{
  e <- environment(lik)
  f.pars <- e$f.pars
  do.asr <- make.do.asr.marginal_divextra(e$all_branches, e$rootfunc)
  asr <- function(pars, nodes = NULL, condition.surv = TRUE,
                  root = ROOT.FLAT, root.p = NULL, ...) {
    pars2 <- f.pars(pars)
    do.asr(pars2, nodes, NULL, condition.surv, root, root.p,
           intermediates = FALSE)
  }
  asr
}



asr.marginal.classe <- function (lik, pars, nodes = NULL, ...) 
  make.asr.marginal.classe(lik)(pars, nodes, ...)


#-----------  ASR along branch

#diversitree:::do.asr.marginal.R
## Next, the utility functions for the different types of models This
## is to asr.marginal what all.branches is for the core models.  Here,
## the argument 'res' is the result of running all.branches
do.asr.marginal.R_branch <- function(pars, cache, res, nodes, states.idx,
                              initial.conditions, branches, root, relative.t.branch, ...) {
  ## Store these for easier calculation.
  children <- cache$children
  parent <- cache$parent
  len <- cache$len
  depth <- cache$depth
  root.idx <- cache$root
  anc <- cache$ancestors

  # Check that nodes length should be one
  if (length(nodes) != 1) {
    stop("Provide a node defining a branch!")
  }
  # nodes <- nodes + cache$n.tip
  nd <- nodes + cache$n.tip
  cat('Working on a branch of the node: ', nd, '\n')

# nd <- 9
#f <- function(nd) {
    cat('##################---- Node: ', nd, '\n')
    ## Include current node but omit root:
    anc.nd <- c(nd, anc[[nd]])
    anc.nd <- anc.nd[-length(anc.nd)]
    p <- rep(NA, length(states.idx))

    # st=1
    for (st in seq_along(states.idx) ) {
      cat('---- st: ', st, '\n')
      lq <- res$lq
      branch.init <- res$init
      branch.base <- res$base
      # !!
      # branch.init[states.idx[-st],nd] <- 0
      y.in <- branch.init[,nd]
      cat('y.in <- branch.init[,nd]: ', y.in, '\n')

      # i=9
      #for ( i in anc.nd ) {
      for (index in 1:length(anc.nd) ) {
        i <- anc.nd[index]

        cat('---- Ancestral Nodes, anc.nd i: ', i, '\n')
        # the first branch is the focal one
        # so ASR calculation is different
        if (index==1){
          # we hame multiple IFs: for numeric and "BASE"
          if (is.numeric(relative.t.branch)){
            # there two ans (from 0 to relative.t.branch) and (from  relative.t.branch to len[i])
            # The Ancestral Probs are at relative.t.branch
            # ans <- branches(y.in, len[i], pars, depth[i], i)
            ans_test <- branches(y.in, len[i], pars, depth[i], i)
            cat('depth[i]: ', depth[i], '\n')
            print('ans_test:')
            print(ans_test)
            # relative.t.branch <- .01
            cat('relative.t.branch: ', relative.t.branch, '\n')
            ans_1 <- branches(y.in, relative.t.branch, pars, depth[i], i)
            #cat('ans_1: ',, '\n')
            print('ans_1:')
            print( ans_1)
            #
            ans_1.lq <- ans_1[[1]]
            ans_1.base <- ans_1[[2]]
            #
            # y.in_2 is where Ancestral Probs are calculated
            y.in_2 <- c(ans_1.base)
            #-----
            y.in_2[states.idx[-st]] <- 0
            cat('y.in_2 <-0: ', y.in_2, '\n')
            ans_2 <- branches(y.in_2, len[i] - relative.t.branch, pars, depth[i] + relative.t.branch, i)
            # ans_2 <- branches(y.in_2, len[i] - relative.t.branch, pars, depth[i] + relative.t.branch, i)
            #cat('ans_2: ', ans_2, '\n')
            print('ans_2:')
            print(ans_2)
            #----
            ans_2.lq <- ans_2[[1]]
            ans_2.base <- ans_2[[2]]
            # we need to sum lqs
            ans_total.lq <- ans_1.lq + ans_2.lq
            # lq[i] <- ans[[1]]
            # branch.base[,i] <- ans[[2]]
            lq[i] <- ans_total.lq
            branch.base[,i] <- ans_2.base
            #
          } else if (relative.t.branch == "BASE") {
            print("BASE")

            branch.base[states.idx[-st],nd] <- 0
            y.in <- branch.base[ ,nd]
            print(y.in)
            lq[i] <- log(sum(y.in))
            print(lq[i])
            #branch.base[,i] <- ans[[2]]

          } else {
            stop("relative.t.branch: arguments has wrong value")
          }

        } else {
          ans <- branches(y.in, len[i], pars, depth[i], i)
          lq[i] <- ans[[1]]
          branch.base[,i] <- ans[[2]]
        }
        #
        j <- parent[i]
        y.in <- initial.conditions(branch.base[,children[j,]], pars,
                                   depth[j], j)
        branch.init[,j] <- y.in
        cat('branch.init[,j] <- y.in: ',  y.in, '\n')
      }
      ans <- root(list(vals=branch.init[,root.idx], lq=lq), pars)
      ## explots IEEE arithmetic's exp(-Inf) == 0
      p[st] <- if ( is.na(ans) ) -Inf else ans
    }

    pp <- exp(p - max(p))
    #pp/sum(pp)
    matrix(pp/sum(pp), ncol = 1)
#}

  #matrix(unlist(lapply(nodes, f)), ncol=length(nodes))
}



make.do.asr.marginal_divextra_branch <- function(all_branches, rootfunc)
{
  eb <- environment(all_branches)
  cache <- eb$cache
  states.idx <- cache$info$idx.d
  if (isTRUE(cache$info$partitioned)) {
    branches <- eb$branches.split
    initial.conditions <- eb$initial.conditions.split
  }
  else {
    # branches <- eb$branches
    branches <- eb$branches.td
    #initial.conditions <- eb$initial.conditions
    initial.conditions <- eb$initial.conditions.td
  }
  function(pars, nodes, preset, relative.t.branch, ...) {
    root.f <- function(res, pars) rootfunc(res, pars, ...)
    res <- all_branches(pars, TRUE, preset)
    do.asr.marginal.R_branch(pars, cache, res, nodes, states.idx,
                             initial.conditions, branches, root.f, relative.t.branch)
  }
}



make.asr.marginal.classe_branch <- function (lik, ...)
{
  e <- environment(lik)
  f.pars <- e$f.pars
  do.asr <- make.do.asr.marginal_divextra_branch(e$all_branches, e$rootfunc)
  asr <- function(pars, nodes, condition.surv = TRUE,
                  root = ROOT.FLAT, root.p = NULL, relative.t.branch,  ...) {
    pars2 <- f.pars(pars)
    do.asr(pars2, nodes, NULL, relative.t.branch, condition.surv, root, root.p,
           intermediates = FALSE)
  }
  asr
}
