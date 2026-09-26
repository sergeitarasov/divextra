#' Remove all simulated extra lineages from a stochastic map
#'
#' @param sim One result from [make.simmap.classe.td()] or
#'   [keep.lineage.simmap()]. For multiple maps use `lapply`.
#' @return A list containing `tree`, `tips`, `lineages`, `events`, `node.map`,
#'   and `selected.tips`, with only the retained observed backbone lineages.
#' @details
#' Removes generated lineages regardless of whether they are extinct or extant.
#' Observed tips are retained. Exact branch maps, ages and the original root
#' are preserved. Unary nodes, including the attachment points of removed
#' side lineages, remain: their state changes at speciation must not be
#' mistaken for anagenetic events. Apply [keep.lineage.simmap()] first to
#' restrict the observed backbone to a focal set of tips.
#' This is pruning of an existing simulated history, not resimulation or
#' conditioning on absence of extra extant taxa.
#' @export
prune.unobserved.simmap <- function(sim) {
  if (is.null(sim$lineages$generated) || is.null(sim$tips$generated) ||
      is.null(sim$tree$maps) || is.null(sim$node.map))
    stop("sim must be one mapped simulation or focused result")
  observed <- sim$tips$label[!sim$tips$generated]
  if (!length(observed)) stop("sim has no observed tips to retain")
  out <- .classe_td_subset_lineages(sim, !sim$lineages$generated, observed)
  # Side-branch transitions are logged at their parent graph node, which may
  # survive pruning. Backbone transitions always carry a backbone.edge ID.
  out$events <- out$events[!(out$events$event == "transition" &
                             is.na(out$events$backbone.edge)), , drop = FALSE]
  out$tree$event.history <- out$events
  out$pruned.unobserved <- TRUE
  out
}
