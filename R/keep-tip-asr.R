#' Keep selected tips and their existing ancestral reconstruction
#'
#' @param x A list returned by [extract.clade.asr()] or `keep.tip.asr()`.
#' @param tips Character vector of exact tip labels to retain (at least two).
#' @return A list with `tree`, `node.probs`, `tip.states`, and `node.map`,
#'   suitable for further pruning or [make.simmap.classe.td()].
#' @details
#' Pruning collapses single-child nodes and sums their branch lengths. Only
#' ASR columns for surviving branching nodes are retained and renumbered;
#' `node.map` continues to refer to the original full-tree node numbers.
#' Probabilities are not recalculated: they still include evidence from taxa
#' removed from the tree. Suppressed nodes and their ASR constraints are not
#' retained along merged branches. Mapping starts at the retained tips' MRCA,
#' without a stem edge. This is not a reconstruction conditional only on the
#' retained tips, nor a mapped history of the removed speciation events.
#' @export keep.tip.asr
keep.tip.asr <- function(x, tips) {
  if (!is.list(x) || !inherits(x$tree, "phylo") ||
      !is.matrix(x$node.probs) || is.null(x$node.map))
    stop("x must be an extract.clade.asr or keep.tip.asr result")
  if (!is.character(tips) || length(tips) < 2L || anyNA(tips) ||
      anyDuplicated(tips))
    stop("tips must contain at least two unique, nonmissing tip labels")
  missing <- setdiff(tips, x$tree$tip.label)
  if (length(missing)) stop("Unknown tip labels: ", paste(missing, collapse = ", "))
  tree <- x$tree
  old.labels <- tree$node.label
  tree$node.label <- as.character(ape::Ntip(tree) + seq_len(tree$Nnode))
  pruned <- ape::keep.tip(tree, tips)
  old.nodes <- as.integer(pruned$node.label)
  new.nodes <- ape::Ntip(pruned) + seq_len(pruned$Nnode)
  if (!all(as.character(old.nodes) %in% colnames(x$node.probs)))
    stop("node.probs is missing retained node columns")
  idx <- match(old.nodes, x$node.map$subtree.node)
  if (anyNA(idx)) stop("node.map is missing retained nodes")
  probs <- x$node.probs[, as.character(old.nodes), drop = FALSE]
  colnames(probs) <- as.character(new.nodes)
  states <- if (is.matrix(x$tip.states)) {
    x$tip.states[, pruned$tip.label, drop = FALSE]
  } else x$tip.states[pruned$tip.label]
  pruned$node.label <- if (is.null(old.labels)) NULL else
    old.labels[old.nodes - ape::Ntip(tree)]
  pruned$root.edge <- NULL
  list(tree = pruned, node.probs = probs, tip.states = states,
       node.map = data.frame(subtree.node = new.nodes,
                             original.node = x$node.map$original.node[idx]))
}
