#' ASR-weighted or tree-conditioned stochastic mapping under episodic ClaSSE
#'
#' @param tree Rooted binary ultrametric backbone, without a stem/root edge.
#' @param pars Named full parameter vector in [make.classe.td()] order:
#'   epoch 1 is nearest the present, with boundaries named `t.1`, `t.2`, etc.
#' @param tip.weights Named tip states or a k-by-Ntip matrix of nonnegative
#'   state probabilities. Supply state information only, without sampling
#'   fractions. The function cannot remove sampling fractions already mixed
#'   into a supplied matrix.
#' @param node.probs Named internal-node states or a k-by-Nnode matrix of
#'   nonnegative ASR weights. Columns use ape node numbers. The root column
#'   supplies the default root distribution in ASR-weighted mode. Omit in
#'   tree-conditioned mode, which computes descendant likelihoods instead.
#' @param root.p Optional nonnegative vector of length k overriding the root
#'   ASR distribution. By default, the root column of `node.probs` is used.
#'   The root state is sampled once per map.
#' @param k Number of states.
#' @param state.labels State labels in parameter order.
#' @param nsim Number of maps.
#' @param max.branch.tries Maximum proposals of a node's daughter pair and
#'   both outgoing branches. Exhaustion raises an error, without resampling
#'   the root or silently discarding a requested map.
#' @param max.events,max.taxa Event and terminal-count limits per map.
#' @param seed Optional random seed.
#' @param side.lineage.tips Optional unique observed tip labels. Expand side
#'   trees only at hidden speciation events on the union of their paths to the
#'   root. NULL (default) expands side trees on every backbone branch. Keep the
#'   complete intended backbone, including any unselected sister tips, as input.
#' @param side.lineage.mode Either `"all"` (default, unrestricted side trees)
#'   or `"extinct_clades"` (each expanded side clade must die out completely).
#' @param max.side.tries Maximum attempts per expanded side clade to achieve
#'   extinction when `extinction.sampler="rejection"`. Exhaustion fails the map;
#'   no side lineage is silently omitted. Unused by the direct sampler.
#' @param extinction.sampler In extinct-clade mode, `"direct"` (default) samples
#'   directly conditional on extinction before present using shared extinction
#'   probability curves. `"rejection"` retains unrestricted forward proposals
#'   followed by rejection of survivors. Ignored in `side.lineage.mode="all"`.
#' @param extinction.control Named list of numerical controls for the direct
#'   extinction-probability solver and event-time inversion; see details.
#' @param backbone.method `"asr-weighted"` preserves the original augmentation
#'   procedure. `"tree-conditioned"` samples backbone histories conditional on
#'   the supplied observed tree and tip state likelihoods, assuming all additional
#'   lineages are extinct at present. This mode defaults to `"extinct_clades"`
#'   when side.lineage.mode is omitted and requires the direct sampler.
#' @param root.prior Explicit nonnegative root-state prior, length k, required
#'   only in tree-conditioned mode. It is multiplied by root descendant
#'   likelihoods. This is NOT a reconstructed root posterior; do not supply ASR
#'   probabilities as a prior unless that is deliberately your model.
#'
#' @details
#' The default implements an ASR-weighted augmentation process, not the joint
#' posterior of ClaSSE conditional on the observed tree. Backbone extinction
#' is disabled by construction. Hidden speciation uses all lambda outcomes;
#' a random daughter continues the backbone. At observed nodes, a daughter
#' pair is drawn proportional to lambda with equal probabilities for its two
#' orientations. Both outgoing branch proposals are discarded if either
#' endpoint fails. Once accepted, these branches are fixed while descendant
#' nodes are processed. Endpoint acceptance uses the supplied weight divided
#' by its maximum. Consequently an accepted endpoint has distribution
#' proportional to its proposal probability times that weight, not generally
#' the supplied ASR marginal itself. No additional Q compatibility factor or
#' root speciation weighting is applied.
#'
#' In `backbone.method="tree-conditioned"`, postorder pruning computes
#' branch-specific descendant likelihoods D alongside the shared E curves.
#' Root states are drawn from root.prior * D(root). Observed daughter states
#' are drawn jointly using lambda and both daughter likelihoods. Along branches,
#' transitions have conditional weights q_ij D_j and the two hidden-speciation
#' orientations have weights lambda_ijk D_j E_k and lambda_ijk D_k E_j.
#' No-event probabilities use exp(-integrated original rate) * D(younger)/D(older),
#' including the nonzero probability of reaching a child without another event.
#' This targets the tree-conditioned extinction-only history distribution up to
#' numerical approximation. Marginal node ASRs and root.p are rejected in this
#' mode to prevent confusing marginal/posterior weights with likelihoods/priors.
#' This is a complete-census target, not the posterior of an incompletely sampled
#' model. Supply the full intended observed tree; deleted observed taxa are not
#' extinct. Use side.lineage.tips and focus pruning after simulation instead.
#'
#' Unrestricted waiting times are exponential and stop at epoch boundaries. After the
#' backbone is accepted, side lineages evolve with the full lambda, mu and Q
#' process. Their initial branches are retained even if they never split.
#' Surviving side tips have fate `extra_extant`; extinct tips have `extinct`.
#' Sampling fractions are not used in endpoint acceptance. They may have
#' informed the externally computed ASR. Side-lineage survivors are retained
#' without conditioning on whether they would have been sampled.
#'
#' With `side.lineage.mode="extinct_clades"`, each expanded side tree is
#' conditioned to have no living descendants at the present. The direct sampler
#' solves the state-dependent extinction equations once per call, from E(0)=0,
#' sharing the resulting curves across every lineage and all `nsim` maps.
#' Epochs retain their original absolute ages and rates. Event times are drawn
#' by inverting exp(-integrated original total rate) * E(younger)/E(older).
#' Event weights are mu, q_ij E_j, and lambda_ijk E_j E_k. This implements the
#' extinction-conditioned process up to numerical solver/interpolation error.
#' Likelihood E values based on incomplete sampling are not used: those describe
#' no sampled descendants rather than literal extinction. The reference rejection
#' sampler retries side trees from the same founding state and age, stopping at
#' the first survivor. Neither sampler forcibly kills survivors or changes
#' observed backbone histories. Hidden births are allowed at all positive ages.
#' In ASR-weighted mode this is only local side-clade extinction conditioning,
#' not a joint ClaSSE posterior: founding ages/states are not reweighted by
#' their extinction probabilities. Tree-conditioned mode also reweights the
#' backbone process using E and D as described above.
#' Rejection may be expensive or impossible. Resource limits fail explicitly,
#' rather than treating oversized proposals as statistical rejections. Event and
#' taxon counters roll back between rejected proposals; limits are not CPU budgets.
#'
#' `extinction.control` accepts `rtol` (1e-10) and `atol` (1e-13) for the
#' extinction ODE, `max.step` (0.01 in tree time units) for the interpolation
#' grid, and `root.tol` (1e-10) for relative event-age inversion accuracy.
#' `interpolation.tol` (1e-6) controls adaptive absolute/log-probability checks;
#' `max.points` (250000) limits each numerical cache grid. Checks of relative
#' probabilities apply above the solver resolution threshold (100 * atol).
#' These are convergence diagnostics, not rigorous global error bounds.
#' Descendant likelihood caches inherit the E solver precision settings.
#' Their log-probability checks allow ODE uncertainty of
#' `10 * rtol + 10 * atol / p` for normalized probability `p`.
#' Per-branch diagnostics identify checks limited by solver precision.
#' Additional geometric grid points resolve ages near zero and epoch boundaries.
#' Linear interpolation adds approximation beyond the ODE tolerance: halve
#' `max.step` and tighten `interpolation.tol` and ODE tolerances to check
#' numerical convergence.
#' Zero extinction probabilities are never replaced with a positive floor.
#' Structurally impossible extinction and inadequate numerical resolution raise
#' errors, rather than silently switching samplers or omitting lineages.
#'
#' Restricting `side.lineage.tips` does not remove observed backbone branches
#' or alter root speciation, daughter orientation, endpoint acceptance, or
#' hidden speciation on unselected branches. Only the subsequent independent
#' side-tree expansion is skipped. All descendants of an expanded side tree
#' are simulated, regardless of their state or eventual fate. Unexpanded side
#' births remain in the event history and may appear as unary nodes. They are
#' not extinctions. Use [keep.lineage.simmap()] afterward to extract the focus.
#' This preserves the focused distribution of the unrestricted simulator when
#' both can complete without resource limits, not necessarily the same
#' realization for a given seed. In a batch, omitted random draws can also
#' change later maps' backbone realizations. This is not additional statistical
#' conditioning and does not make the ASR-weighted sampler an exact posterior.
#'
#' @return A `classe_td_tree_sim` with `tree` (a simmap), `tips`, `events`,
#'   `lineages`, `node.map`, `branch.attempts`, and method metadata. For nsim
#'   greater than one, a `classe_td_tree_simulations` list. Plot with
#'   `phytools::plotSimmap(sim$tree, colors=colors, split.vertical=TRUE)`.
#'   `side.lineage.edges` records eligible input-tree edge rows;
#'   `side.lineage.events` records every accepted hidden birth and whether its
#'   side tree was expanded. Counts outside selected paths are not complete.
#'   `side.attempts` gives attempts for each expanded hidden-birth graph node.
#'   `side.lineage.mode`, `extinction.sampler`, and `extinction.control` record
#'   the conditioning and numerical settings. `side.attempts` is one per expanded
#'   side clade for the direct sampler.
#'   In tree-conditioned mode, root.prior and root.loglik record the prior and
#'   unconditioned-on-crown-survival log likelihood (using normalized tip weights).
#'   `conditioned=TRUE` and `numerical.approximation=TRUE` distinguish this target
#'   from the original augmentation. `posterior.exact` remains FALSE because the
#'   implementation uses numerical integration/interpolation.
#' @export
make.simmap.classe.td <- function(
    tree, pars, tip.weights, node.probs = NULL, root.p = NULL, k,
    state.labels = as.character(seq_len(k)), nsim = 1L,
    max.branch.tries = 10000L, max.events = 100000L,
    max.taxa = 100000L, seed = NULL, side.lineage.tips = NULL,
    side.lineage.mode = c("all", "extinct_clades"), max.side.tries = 10000L,
    extinction.sampler = c("direct", "rejection"), extinction.control = list(),
    backbone.method = c("asr-weighted", "tree-conditioned"), root.prior = NULL) {
  backbone.method <- match.arg(backbone.method)
  tree.conditioned <- backbone.method == "tree-conditioned"
  if (tree.conditioned && missing(side.lineage.mode)) side.lineage.mode <- "extinct_clades"
  side.lineage.mode <- match.arg(side.lineage.mode)
  extinction.sampler <- match.arg(extinction.sampler)
  if (tree.conditioned) {
    if (!is.null(node.probs) || !is.null(root.p))
      stop("tree-conditioned mode uses descendant likelihoods, not node.probs/root.p; supply root.prior")
    if (side.lineage.mode != "extinct_clades" || extinction.sampler != "direct")
      stop("tree-conditioned mode requires extinct_clades with the direct sampler")
    if (is.null(root.prior)) stop("Supply an explicit root.prior in tree-conditioned mode")
    node.probs <- matrix(1, k, tree$Nnode, dimnames=list(state.labels,
      as.character(length(tree$tip.label)+seq_len(tree$Nnode))))
    root.p <- root.prior
  } else {
    if (is.null(node.probs)) stop("Supply node.probs in asr-weighted mode")
    if (!is.null(root.prior)) stop("root.prior is only used in tree-conditioned mode")
  }
  positive.integer <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) &&
      x >= 1 && x == floor(x)
  }
  for (nm in c("k", "nsim", "max.branch.tries", "max.events", "max.taxa", "max.side.tries")) {
    if (!positive.integer(get(nm))) stop(nm, " must be a positive integer")
  }
  if (length(state.labels) != k || anyNA(state.labels) ||
      anyDuplicated(state.labels) || any(!nzchar(state.labels)))
    stop("state.labels must contain k unique non-empty labels")
  if (!is.null(tree$root.edge) && tree$root.edge != 0)
    stop("This version starts at the crown root; supply a tree without root.edge")
  if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) ||
      any(tree$edge.length <= 0))
    stop("Backbone edge lengths must be finite and strictly positive")
  backbone <- mark.classe.td.backbone(
    tree, node.probs, tip.weights, k, state.labels
  )
  side.edges <- seq_len(nrow(tree$edge))
  if (!is.null(side.lineage.tips)) {
    if (!is.character(side.lineage.tips) || !length(side.lineage.tips) ||
        anyNA(side.lineage.tips) || anyDuplicated(side.lineage.tips) ||
        !all(side.lineage.tips %in% tree$tip.label))
      stop("side.lineage.tips must contain unique observed tip labels")
    eligible <- rep(FALSE, nrow(tree$edge))
    for (node in match(side.lineage.tips, tree$tip.label)) {
      while (node != backbone$root) {
        row <- match(node, tree$edge[, 2])
        eligible[row] <- TRUE
        node <- tree$edge[row, 1]
      }
    }
    side.edges <- which(eligible)
  }
  root.source <- if (is.null(root.p)) "node-ASR" else "custom"
  if (is.null(root.p))
    root.p <- backbone$node.probs[, as.character(backbone$root)]
  if (!is.numeric(root.p) || length(root.p) != k ||
      any(!is.finite(root.p)) || any(root.p < 0) || sum(root.p) <= 0)
    stop("root.p must contain k finite nonnegative weights with positive sum")
  if (!is.null(names(root.p))) {
    if (!setequal(names(root.p), state.labels))
      stop("Named root.p must match state.labels")
    root.p <- root.p[state.labels]
  }
  root.p <- root.p / sum(root.p)
  # Columnwise constants cancel from accepted distributions. Max scaling is
  # faster than sum scaling and retains relative state weights.
  backbone$node.probs <- sweep(backbone$node.probs, 2L,
                               apply(backbone$node.probs, 2L, max), "/")
  backbone$tip.probs <- sweep(backbone$tip.probs, 2L,
                              apply(backbone$tip.probs, 2L, max), "/")
  schedule <- .classe_td_schedule(pars, k, state.labels)
  if (any(schedule$boundaries < 0)) stop("Epoch boundaries must be nonnegative")
  extinction.cache <- if (side.lineage.mode == "extinct_clades" &&
                           extinction.sampler == "direct")
    .classe_td_extinction_cache(schedule, backbone$height, extinction.control) else NULL
  descendant.cache <- if (tree.conditioned)
    .classe_td_descendant_cache(backbone, schedule, extinction.cache) else NULL
  if (tree.conditioned) {
    root.prior <- root.p
    root.p <- root.prior * descendant.cache$root
    if (sum(root.p) <= 0) stop("Observed tree has zero likelihood under the supplied root prior")
    root.loglik <- log(sum(root.p)) + descendant.cache$root.log.scale
    root.p <- root.p / sum(root.p)
    root.source <- "prior-times-descendant-likelihood"
  }
  if (!is.null(seed)) set.seed(seed)
  simulations <- vector("list", nsim)
  for (i in seq_len(nsim)) {
    root.state <- sample.int(k, 1L, prob = root.p)
    ans <- .simulate_classe_td_once(
      backbone, schedule, max.events, max.taxa, max.branch.tries,
      guided = TRUE, initial.root.state = root.state,
      side.lineage.edges = if (is.null(side.lineage.tips)) NULL else side.edges,
      side.lineage.mode = side.lineage.mode, max.side.tries = max.side.tries,
      extinction.cache = extinction.cache, descendant.cache = descendant.cache
    )
    if (!isTRUE(ans$success))
      stop("Simulation ", i, " failed: ", ans$reason,
           ". No completed map is returned for this call.")
    ans$backbone <- backbone
    ans$parameters <- pars
    ans$root.p <- root.p
    ans$root.source <- root.source
    ans$method <- "asr-weighted-backbone"
    ans$posterior.exact <- FALSE
    ans$conditioned <- FALSE
    ans$endpoint.constraint <- "relative-weights"
    if (tree.conditioned) {
      ans$root.prior <- root.prior
      ans$root.loglik <- root.loglik
      ans$method <- "tree-conditioned-extinction-only"
      ans$conditioned <- TRUE
      ans$endpoint.constraint <- "descendant-likelihood"
      ans$conditioning.target <- "supplied tree and tip states; no additional survivors"
      ans$numerical.approximation <- TRUE
      ans$numerical.diagnostics <- list(
        extinction=extinction.cache$diagnostics,
        descendant.points=vapply(descendant.cache$branches,
          function(x)length(x$times),integer(1)),
        descendant=lapply(descendant.cache$branches,function(x)x$diagnostics))
    }
    ans$simulation <- i
    ans$side.lineage.tips <- side.lineage.tips
    ans$side.lineage.edges <- side.edges
    ans$side.lineage.mode <- side.lineage.mode
    if (side.lineage.mode == "extinct_clades") {
      ans$extinction.sampler <- extinction.sampler
      if (!is.null(extinction.cache)) ans$extinction.control <- extinction.cache$control
    }
    births <- ans$events[ans$events$event == "hidden_speciation", , drop = FALSE]
    ans$side.lineage.events <- data.frame(
      event.id = births$event.id, node = births$node, age = births$age,
      backbone.edge = births$backbone.edge, state = births$to2.label,
      expanded = births$backbone.edge %in% side.edges
    )
    class(ans) <- "classe_td_tree_sim"
    simulations[[i]] <- ans
  }
  .classe_td_simulation_set(simulations)
}
