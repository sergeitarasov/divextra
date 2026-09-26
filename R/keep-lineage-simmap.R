#' Retain selected backbone paths and their simulated side lineages
#'
#' @param sim One result from [make.simmap.classe.td()], or a batch returned
#'   with `nsim > 1`. A nonempty ordinary list of simulations is also accepted.
#' @param tips One or more exact observed backbone tip labels.
#' @return A list containing `tree` (a stochastic mapping tree), `tips`,
#'   `lineages`, `events`, `node.map`, and `selected.tips` for a single input.
#'   For a batch, returns a list of these results in input order, preserving
#'   list names. A one-element input list remains a one-element output list.
#' @details
#' Keeps the union of paths from selected observed tips to the original root,
#' plus all generated subtrees born on those paths. Other observed branches,
#' including the other root daughter, and their generated subtrees are removed.
#' Unary nodes are deliberately retained to preserve ages, mapped segments and
#' the original root. Events at retained branching nodes describe the original
#' simulation; a daughter may have been pruned. This is a filtered history,
#' not a new simulation conditional on the selected tips.
#' @export
keep.lineage.simmap <- function(sim, tips) {
  is.single <- function(x) is.list(x) && inherits(x$tree, "phylo") &&
    !is.null(x$lineages) && !is.null(x$node.map) && !is.null(x$tree$maps)
  if (!is.single(sim)) {
    if (!is.list(sim) || !length(sim))
      stop("sim must be one mapped simulation or a nonempty list of simulations")
    valid <- vapply(sim, is.single, logical(1))
    if (!all(valid)) stop("Invalid mapped simulation at batch index ", which(!valid)[1])
    out <- lapply(seq_along(sim), function(i) {
      tryCatch(keep.lineage.simmap(sim[[i]], tips), error = function(e)
        stop("Simulation ", i, ": ", conditionMessage(e), call. = FALSE))
    })
    names(out) <- names(sim)
    return(out)
  }
  observed <- sim$tips$label[!sim$tips$generated]
  if (!is.character(tips) || !length(tips) || anyNA(tips) || anyDuplicated(tips) ||
      !all(tips %in% observed))
    stop("tips must be unique observed backbone tip labels")
  lin <- sim$lineages
  keep <- rep(FALSE, nrow(lin))
  for (node in sim$tips$node[match(tips, sim$tips$label)]) {
    repeat {
      i <- match(node, lin$child)
      if (is.na(i)) break
      keep[i] <- TRUE
      node <- lin$parent[i]
    }
  }
  # Follow only generated edges off the selected paths, never other backbone
  # daughters of an observed node (in particular, the other root daughter).
  repeat {
    nodes <- unique(c(lin$parent[keep], lin$child[keep]))
    more <- lin$generated & lin$parent %in% nodes
    if (all(keep | !more)) break
    keep <- keep | more
  }
  .classe_td_subset_lineages(sim, keep, tips)
}

.classe_td_subset_lineages <- function(sim, keep, tips = NULL) {
  lin <- sim$lineages[keep, , drop = FALSE]
  tip.table <- sim$tips[sim$tips$node %in% lin$child, , drop = FALSE]
  tip.list <- lapply(seq_len(nrow(tip.table)), function(i) tip.table[i, , drop = FALSE])
  names(tip.list) <- tip.table$node
  ids <- setNames(sim$node.map$phylo.node, sim$node.map$graph.node)
  edge.keys <- paste(sim$tree$edge[, 1], sim$tree$edge[, 2])
  edges <- lapply(seq_len(nrow(lin)), function(i) {
    x <- lin[i, ]
    j <- match(paste(ids[x$parent], ids[x$child]), edge.keys)
    map <- if (isTRUE(x$root.edge)) sim$tree$root.map else sim$tree$maps[[j]]
    list(parent = x$parent, child = x$child, start.age = x$start.age,
         end.age = x$end.age, length = x$length, map = map,
         backbone.edge = x$backbone.edge, generated = x$generated)
  })
  out <- .classe_td_graph_to_simmap(edges, tip.list, colnames(sim$tree$mapped.edge),
    root.start = if (any(lin$root.edge)) lin$parent[which(lin$root.edge)[1]] else NULL)
  nodes <- unique(c(lin$parent, lin$child))
  events <- sim$events[sim$events$node %in% nodes &
    (is.na(sim$events$backbone.edge) |
       sim$events$backbone.edge %in% lin$backbone.edge), , drop = FALSE]
  events$phylo.node <- out$node.map$phylo.node[match(events$node, out$node.map$graph.node)]
  out$events <- events
  out$selected.tips <- tips
  out$tree$tip.status <- out$tips
  out$tree$node.map <- out$node.map
  out$tree$lineage.history <- out$lineages
  out$tree$event.history <- events
  out$tree$cladogenetic.events <- events[events$event %in%
    c("observed_speciation", "hidden_speciation", "speciation"), , drop = FALSE]
  out
}
