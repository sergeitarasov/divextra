

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


#' @export
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


#' Ancestral State Reconstruction for ClaSSE/GeoSSE Models
#'
#' Perform marginal ancestral state reconstruction for classe models. This function
#' computes the marginal probability of each state at ancestral nodes given the
#' observed tip states and model parameters.
#'
#' @param lik A likelihood function created by \code{make.classe} or related functions.
#' @param pars A vector of parameters suitable for the likelihood function.
#' @param nodes An optional vector of nodes to return ancestral states for (using ape's
#'   node indexing). By default, all internal nodes are returned.
#' @param ... Additional arguments passed through to the reconstruction function,
#'   including root conditions and survival conditioning.
#'
#' @details
#' This function performs marginal ancestral state reconstruction, computing the
#' probability distribution over character states at each ancestral node. The
#' reconstruction takes into account the phylogenetic relationships, branch lengths,
#' and the evolutionary model specified in the likelihood function.
#'
#' The function is a wrapper that creates and calls the reconstruction function
#' using \code{make.asr.marginal.classe}.
#'
#' @return A matrix where each column represents a node and each row represents
#'   a character state. The values are the marginal probabilities of each state
#'   at each node.
#'
#' @seealso \code{\link{make.asr.marginal.classe}} for the function factory version,
#'   \code{asr.marginal} in the diversitree package for the generic function.
#'
#' @examples
#' #------- Toy dataset, parameters given
#' data(phy5)
#' cols <- plyr::mapvalues(phy5$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
#' plot(phy5, label.offset = 2, cex=1, no.margin = FALSE)
#' tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
#' nodelabels()
#' tiplabels()
#' edgelabels()
#' axisPhylo()
#'
#' lik.td <-make.classe.td(phy5, phy5$tip.state_3, k=3, n.epoch=2, control=list(backend = "gslode"),
#'  strict=TRUE)
#' par.td <- c(4, c(1:54)/100)
#' names(par.td) <- diversitree::argnames(lik.td)
#' lik.est <- lik.td(par.td, intermediates=TRUE)
#' print(lik.est)
#'
#' st <-  asr.marginal.classe(lik.td, par.td)
#' nodelabels(thermo=t(st), piecol=c('black', "blue","red"), cex=2, adj=.5)
#'
#'
#------- Toy dataset, parameters estimated
#' \dontrun{
#' file_path <- system.file("extdata", "geosse_tb_tree.rds", package = "divextra")
#' phy <- readRDS(file_path)
#' plot(phy)
#'
#' file_path <- system.file("extdata", "geosse3_td_test.yml", package = "divextra")
#' par.categories.td <- read_yaml_pars_td(file_path)
#' # display
#' print(par.categories.td)
#'
#' lik.td <-make.classe.td(phy, phy$tip.state, k=3, n.epoch=2, control=list(backend = "gslode"),
#'  strict=TRUE)
#' formula.td <- make_constraints_sse_td(par.categories.td)
#' lik.const.td <- constrain(lik.td, formulae = formula.td)
#' starting.point <- init.pars.classe_td(lik.const.td, phy, k=3, n.epoch=2, eps=0.5)
#' print(starting.point)
#'
#' mle.td <- find.mle(lik.const.td, starting.point, condition.surv=TRUE, keep.func=FALSE)
#'
#' st <- asr.marginal.classe(lik.td, mle.td$par.full)
#'
#' cols <- plyr::mapvalues(phy$tip.state, from = c(1:3), to=c('orange', "green","red" ))
#' plot(phy, show.tip.label = FALSE, cex=1, no.margin = FALSE)
#' tiplabels(pch = 15, col = cols, cex = .5, adj = 1)
#' nodelabels(thermo=t(st), piecol=c('orange', "green","red" ), cex=.3, adj=0)
#' }
#'
#' @export
asr.marginal.classe <- function (lik, pars, nodes = NULL, ...)
  make.asr.marginal.classe(lik)(pars, nodes, ...)



#-----------------------------------------
#-----------  ASR along branch
#-----------------------------------------

# diversitree:::do.asr.marginal.R
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

  ## Include current node but omit root:
  anc.nd <- c(nd, anc[[nd]])
  anc.nd <- anc.nd[-length(anc.nd)]
  p <- rep(NA, length(states.idx))

  # state st=1
  for (st in seq_along(states.idx) ) {
    lq <- res$lq
    branch.init <- res$init
    branch.base <- res$base
    # branch.init[states.idx[-st],nd] <- 0
    y.in <- branch.init[,nd]

    # current branch index: index=9
    for (index in 1:length(anc.nd) ) {
      i <- anc.nd[index]
      cat('Node: ', i, '\n')
      # the first branch is the focal one -- so ASR calculation is different
      if (index==1){
        # we hame multiple IFs: for numeric and "BASE" arguments of relative.t.branch
        if (is.numeric(relative.t.branch)){
          # there two ans (from 0 to relative.t.branch) and (from  relative.t.branch to len[i])
          # The Ancestral Probs are at relative.t.branch
          ans_1 <- branches(y.in, relative.t.branch, pars, depth[i], i)
          ans_1.lq <- ans_1[[1]]
          ans_1.base <- ans_1[[2]]
          #
          # y.in_2 is where Ancestral Probs are calculated
          y.in_2 <- c(ans_1.base)
          y.in_2[states.idx[-st]] <- 0
          ans_2 <- branches(y.in_2, len[i] - relative.t.branch, pars, depth[i] + relative.t.branch, i)
          ans_2.lq <- ans_2[[1]]
          ans_2.base <- ans_2[[2]]
          # we need to sum lqs
          ans_total.lq <- ans_1.lq + ans_2.lq
          lq[i] <- ans_total.lq
          branch.base[,i] <- ans_2.base

        } else if (relative.t.branch == "BASE") {
          # print("BASE")
          branch.base[states.idx[-st],nd] <- 0
          y.in <- branch.base[ ,nd]
          lq[i] <- log(sum(y.in))
        } else {
          stop("relative.t.branch: arguments has wrong value")
        }
      # if "index" is not ==1
      } else {
        ans <- branches(y.in, len[i], pars, depth[i], i)
        lq[i] <- ans[[1]]
        branch.base[,i] <- ans[[2]]
      }
      #
      j <- parent[i]
      y.in <- initial.conditions(branch.base[,children[j,]], pars, depth[j], j)
      branch.init[,j] <- y.in
    }
    ans <- root(list(vals=branch.init[,root.idx], lq=lq), pars)
    ## explots IEEE arithmetic's exp(-Inf) == 0
    p[st] <- if ( is.na(ans) ) -Inf else ans
  }

  pp <- exp(p - max(p))
  #pp/sum(pp)
  matrix(pp/sum(pp), ncol = 1)
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



#' @export
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




#' Ancestral State Reconstruction Along a Branch for ClaSSE/GeoSSE Models
#'
#' Perform marginal ancestral state reconstruction at a specific time point along
#' a branch for classe models. This function computes the marginal probability of
#' each state at a given relative time along a specific branch.
#'
#' @param lik A likelihood function created by \code{make.classe.td} or related functions.
#' @param pars A vector of parameters suitable for the likelihood function.
#' @param node.id The node ID defining the branch of interest (using ape's node indexing,
#'   same as displayed by \code{ape::nodelabels()}).
#' @param relative.t.branch The relative time along the branch where reconstruction
#'   should be performed. Time starts from present time (0) and increases as
#'   the branch continues. Can be numeric or "BASE" for the base of the branch (i.e.,
#'   where the focal branch joins with its sister branch).
#'   Must be <= branch length.
#' @param ... Additional arguments passed through to the reconstruction function,
#'   including root conditions and survival conditioning.
#'
#' @details
#' This function performs marginal ancestral state reconstruction at a specific
#' point along a branch, allowing examination of evolutionary transitions within
#' branches rather than just at nodes. The reconstruction accounts for the
#' phylogenetic relationships, branch lengths, and time-dependent evolutionary
#' model specified in the likelihood function.
#'
#' The function internally converts node IDs from ape's indexing system to
#' diversitree's indexing system and validates that the requested time point
#' does not exceed the branch length.
#'
#'
#' @return A column vector containing the marginal probabilities of each character
#'   state at the specified time point along the branch.
#'
#' @seealso \code{\link{asr.marginal.classe}} for reconstruction at nodes,
#'   \code{\link{make.asr.marginal.classe_branch}} for the function factory version.
#'
#' @examples
#' data(phy5)
#' cols <- plyr::mapvalues(phy5$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
#' plot(phy5, label.offset = 2, cex=1, no.margin = FALSE)
#' tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
#' nodelabels()
#' tiplabels()
#' edgelabels()
#' axisPhylo()
#'
#' lik.td <-make.classe.td(phy5, phy5$tip.state_3, k=3, n.epoch=2, control=list(backend = "gslode"),
#'   strict=TRUE)
#' par.td <- c(4, c(1:54)/100)
#' names(par.td) <- diversitree::argnames(lik.td)
#' lik.est <- lik.td(par.td, intermediates=TRUE)
#' print(lik.est)
#'
#' st <-  asr.marginal.classe(lik.td, par.td)
#' # ASR reconstruction at t=1 along branch 3 (defined by terminal node 9)
#' # Note: relative.t.branch starts from present time and equals 0 at the branch's start
#' st9 <-asr.marginal.classe_branch(lik.td, par.td, node.id=9, relative.t.branch= 1 )
#'
#' print(st)
#' print(st9)
#'
#' @export
asr.marginal.classe_branch <- function (lik, pars, node.id, relative.t.branch, ...){
  e <- environment(lik)
  eb <- environment(e$all_branches)
  cache <- eb$cache

  # diversitree uses weird node ids
  # convert our nod id to diversitree
  # our node is the same as plotiing with ape::nodelabels()
  node.id.diversitree <- node.id - Ntip(cache$info$phy)

  br.len <- cache$len[node.id]
  # Check that relative.t.branch should be the same or shorter than branch length
  if (is.numeric(relative.t.branch) && (relative.t.branch > br.len)) {
    stop(paste("relative.t.branch (", relative.t.branch,
               ") should be the same or shorter than branch length (", br.len, ")"))
  }

  make.asr.marginal.classe_branch(lik)(pars, nodes=node.id.diversitree, relative.t.branch=relative.t.branch, ...)
}




#' Ancestral State Reconstruction Along a Branch at Multiple Time Points
#'
#' Perform marginal ancestral state reconstruction at multiple time points along
#' a specific branch for classe models. This function computes the marginal
#' probability of each state at evenly spaced time intervals along a branch,
#' providing a detailed view of evolutionary transitions within the branch.
#'
#' @param lik A likelihood function created by \code{make.classe.td} or related functions.
#' @param pars A vector of parameters suitable for the likelihood function.
#' @param node.id The node ID defining the branch of interest (using ape's node indexing,
#'   same as displayed by \code{ape::nodelabels()}).
#' @param Nbins The number of time bins to sample along the branch (default: 1).
#'   More bins provide higher resolution but require more computation.
#' @param eps A small epsilon value to avoid sampling exactly at the branch endpoints
#'   (default: 0.01). This prevents numerical issues at boundary conditions.
#' @param include.ends Logical indicating whether to include reconstructions at the
#'   exact branch endpoints (present time and base of branch). If TRUE, adds
#'   reconstructions at t=0 and t="BASE" to the results (default: TRUE).
#' @param ... Additional arguments passed through to the reconstruction function,
#'   including root conditions and survival conditioning.
#'
#' @details
#' This function extends \code{asr.marginal.classe_branch} to perform reconstructions
#' at multiple time points along a single branch. Time points are evenly distributed
#' between \code{eps} and \code{branch_length - eps}. When \code{include.ends=TRUE},
#' additional reconstructions are added at the exact endpoints.
#'
#' The function is particularly useful for:
#' \itemize{
#'   \item Visualizing evolutionary transitions along branches
#'   \item Detecting rapid evolutionary changes within branches
#'   \item Creating detailed plots of state probability changes over time
#' }
#'
#' Time measurement follows the same convention as \code{asr.marginal.classe_branch}:
#' t=0 corresponds to present time (tip), increasing values move toward the ancestral node.
#'
#' @return A list containing:
#' \describe{
#'   \item{asr}{A matrix where each column represents a time point and each row
#'     represents a character state. Values are marginal probabilities.}
#'   \item{t}{A vector of time points corresponding to the columns in \code{asr}.}
#'   \item{depth}{The depth (distance from root) of the focal node.}
#' }
#'
#' @seealso \code{\link{asr.marginal.classe_branch}} for single time point reconstruction,
#'   \code{\link{asr.marginal.classe}} for reconstruction at nodes only.
#'
#' @examples
#' data(phy5)
#' cols <- plyr::mapvalues(phy5$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
#' plot(phy5, label.offset = 2, cex=1, no.margin = FALSE)
#' tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
#' nodelabels()
#' tiplabels()
#' edgelabels()
#' axisPhylo()
#'
#' lik.td <-make.classe.td(phy5, phy5$tip.state_3, k=3, n.epoch=2, control=list(backend = "gslode"),
#'   strict=TRUE)
#' par.td <- c(4, c(1:54)/100)
#' names(par.td) <- diversitree::argnames(lik.td)
#' lik.est <- lik.td(par.td, intermediates=TRUE)
#' print(lik.est)
#'
#' st <-  asr.marginal.classe(lik.td, par.td)
#' # ASR reconstruction at t=1 along branch 3 (defined by terminal node 9)
#' # Note: relative.t.branch starts from present time and equals 0 at the branch's start
#' br9 <- asr.marginal.classe_branch.multiple(lik.td, par.td, node.id=9, Nbins=50, include.ends=TRUE)
#' plot(br9$t, br9$asr[1,], type='l', xlim=rev(range(br9$t)))
#'
#' @export
asr.marginal.classe_branch.multiple <- function (lik, pars, node.id, Nbins=1, eps=0.01, include.ends=TRUE,  ...){
  e <- environment(lik)
  eb <- environment(e$all_branches)
  cache <- eb$cache


  # diversitree uses weird node ids
  # convert our nod id to diversitree
  # our node is the same as plotiing with ape::nodelabels()
  node.id.diversitree <- node.id - Ntip(cache$info$phy)

  br.len <- cache$len[node.id]
  time.bins <- seq(0+eps, br.len-eps, length.out=Nbins)
  ans.t <- c()
  for (t in time.bins){
    #print(t)
    x <- make.asr.marginal.classe_branch(lik)(pars, nodes=node.id.diversitree, relative.t.branch= t, ...)
    ans.t <- cbind(ans.t, x)
  }


  if (include.ends==TRUE){
    br.0 <- make.asr.marginal.classe_branch(lik)(pars, nodes=node.id.diversitree, relative.t.branch= 0, ...)
    br.1 <- make.asr.marginal.classe_branch(lik)(pars, nodes=node.id.diversitree, relative.t.branch= "BASE", ...)
    ans.t <- cbind(br.0, ans.t, br.1)
    time.bins <- c(0, time.bins, br.len)
  }

  return(list(asr=ans.t, t=time.bins, depth=cache$depth[node.id]))

}
