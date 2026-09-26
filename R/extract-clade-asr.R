#' Extract a clade together with matching ancestral and tip states
#'
#' @param tree Full rooted phylogeny used for ancestral reconstruction.
#' @param node Focal internal node in the full tree (ape numbering).
#' @param asr State-by-node matrix from [asr.marginal.classe()]. Unnamed
#'   columns must include all internal nodes in ape order. Named columns
#'   must use original ape node numbers.
#' @param tip.states Named full-tree tip-state vector, or state-by-tip matrix
#'   with columns named by tip labels. Supply states without sampling fractions.
#' @param state.labels Optional unique labels in the current ASR row order.
#'   If supplied, assigns ASR row names without reordering rows. Otherwise
#'   existing ASR row names are preserved (including NULL if unnamed).
#' @return A list with `tree`, `node.probs`, `tip.states`, and `node.map`.
#'   ASR columns are renamed to the extracted tree's node numbers; `node.map`
#'   records their original numbers. The extracted tree has no stem edge.
#' @description
#' Selects already calculated full-tree ASR without recalculating it on the
#' subclade. Internal nodes are matched using their original node identities.
#' @export
extract.clade.asr <- function(tree, node, asr, tip.states, state.labels = NULL) {
  internal <- ape::Ntip(tree) + seq_len(tree$Nnode)
  if (length(node) != 1L || is.na(node) || !node %in% internal)
    stop("node must identify an internal node of tree")
  if (!is.matrix(asr)) stop("asr must be a state-by-node matrix")
  if (!is.null(state.labels)) {
    if (!is.character(state.labels) || length(state.labels) != nrow(asr) ||
        anyNA(state.labels) || any(!nzchar(state.labels)) || anyDuplicated(state.labels))
      stop("state.labels must contain one unique nonempty label per ASR row")
    rownames(asr) <- state.labels
  }
  if (is.null(colnames(asr))) {
    if (ncol(asr) != tree$Nnode)
      stop("Unnamed asr must contain every internal node")
    colnames(asr) <- as.character(internal)
  }
  if (anyDuplicated(colnames(asr))) stop("ASR node names must be unique")
  original.labels <- tree$node.label
  tree$node.label <- as.character(internal)
  sub <- ape::extract.clade(tree, node)
  original.nodes <- as.integer(sub$node.label)
  sub.nodes <- ape::Ntip(sub) + seq_len(sub$Nnode)
  if (!all(as.character(original.nodes) %in% colnames(asr)))
    stop("asr is missing nodes required by the extracted clade")
  probs <- asr[, as.character(original.nodes), drop = FALSE]
  colnames(probs) <- as.character(sub.nodes)
  if (is.matrix(tip.states)) {
    if (!all(sub$tip.label %in% colnames(tip.states)))
      stop("tip.states is missing clade tip labels")
    tips <- tip.states[, sub$tip.label, drop = FALSE]
  } else {
    if (!all(sub$tip.label %in% names(tip.states)))
      stop("tip.states must name every clade tip")
    tips <- tip.states[sub$tip.label]
  }
  sub$node.label <- if (is.null(original.labels)) NULL else
    original.labels[original.nodes - ape::Ntip(tree)]
  sub$root.edge <- NULL
  list(tree = sub, node.probs = probs, tip.states = tips,
       node.map = data.frame(subtree.node = sub.nodes, original.node = original.nodes))
}
