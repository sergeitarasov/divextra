# Entire 22-tip Sisyphini example from node 495 in sim-on-tree/asr-new.R.
# Run this script from the DivExtra source checkout, where sim-on-tree/ contains
# the fitted model and source tree.

# This example uses functions from the current source checkout. Calling only
# library(divextra) may load an older installed copy of the package.
if (file.exists("R/sim-on-tree-classe-td.R")) {
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Install pkgload or install the current DivExtra checkout first")
  pkgload::load_all(".", export_all = FALSE)
} else {
  library(divextra)
}

phy <- readRDS("sim-on-tree/phy-271.rds")
geo.regions.abc <- readRDS("sim-on-tree/geo_regions-271.rds")
fit <- readRDS("sim-on-tree/AfrSunken16_gr-6_r-3.rds")
n.tip <- ape::Ntip(phy)

descendants <- function(tree, node) {
  out <- integer()
  queue <- node
  while (length(queue)) {
    children <- tree$edge[tree$edge[, 1] %in% queue, 2]
    out <- c(out, children)
    queue <- children[children > ape::Ntip(tree)]
  }
  out
}

# Original node 495 is the entire focal Sisyphini clade.
focal.node <- 495
all.descendants <- descendants(phy, focal.node)
tips <- all.descendants
tips <- tips[tips <= n.tip]
sisyphini <- ape::drop.tip(phy, phy$tip.label[-tips])
plot(sisyphini)

state.labels <- c(
  "A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
  "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R"
)

# Use `st` already created by sim-on-tree/asr-new.R when available. Otherwise,
# reconstruct the full-tree marginal states from the fitted model.
region.codes <- c(Afr = 1L, OP = 2L, Mada = 3L, Aus = 4L, Maur = 6L)
full.tip.states <- unname(region.codes[geo.regions.abc])
names(full.tip.states) <- names(geo.regions.abc)

if (!exists("st", inherits = TRUE)) {
  sampling.f <- c(0.039, 0.027, 0.446, 0.008, 1e-6, 0.833, rep(1e-6, 10))
  names(sampling.f) <- seq_len(16)
  likelihood <- make.classe.td(
    phy, full.tip.states, k = 16, n.epoch = 2,
    control = list(backend = "gslode"), strict = FALSE,
    sampling.f = sampling.f
  )
  st <- asr.marginal.classe(
    likelihood, fit$par.full,
    root = diversitree::ROOT.GIVEN,
    root.p = c(1, rep(0, 15))
  )
}

original.nodes <- sort(c(
  focal.node, all.descendants[all.descendants > n.tip]
))
# Keep the complete marginal reconstruction. At an internal node, a simulated
# arrival in state i is accepted with this node's reconstructed probability for
# state i. This avoids conditioning every branch on a single MAP endpoint.
node.states <- st[, original.nodes - n.tip, drop = FALSE]
colnames(node.states) <- ape::Ntip(sisyphini) + seq_len(sisyphini$Nnode)
tip.states <- full.tip.states[sisyphini$tip.label]

sim <- simulate.extinct.classe.td(
  sisyphini,
  fit$par.full,
  node.states = node.states,
  tip.states = tip.states,
  k = 16,
  state.labels = state.labels,
  max.tries = 100000,
  max.branch.tries = 100,
  # This is an additional visualization condition. Without it, a valid draw
  # may contain no extinct species, especially for a short focal clade.
  min.extinct = 0,
  nsim = 1,
  seed = 10
)

stopifnot(reconstruct.classe.td(sim)$matches.backbone)

cols <- stats::setNames(
  grDevices::hcl.colors(length(state.labels), "Dynamic"), state.labels
)
phytools::plotSimmap(sim$tree, colors = cols, ftype = "i")
legend("topleft", legend = state.labels, fill = cols, cex = 0.6, bty = "n")

sim$tree$tip.status
sim$tree$cladogenetic.events
sim$tree$event.history
