devtools::load_all("/Users/taravser/GitHub/divextra")
library(ape)
library(phytools)

# Run from the DivExtra package root, as with asr-new.R.
phy <- readRDS("sim-on-tree/phy-271.rds")
geo_regions_abc <- readRDS("sim-on-tree/geo_regions-271.rds")
mle.td <- readRDS("sim-on-tree/AfrSunken16_gr-6_r-3.rds")

geo_regions <- setNames(
  unname(c(Afr = 1L, OP = 2L, Mada = 3L, Aus = 4L, Maur = 6L)[
    as.character(geo_regions_abc)]), names(geo_regions_abc)
)
state.labels <- c("A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
                  "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R")




# First calculate ASR for every internal node of the entire tree.
# Sampling fractions enter this likelihood only, not mapping acceptance.
lik.td <- make.classe.td(
  phy, geo_regions, k = 16, n.epoch = 2, strict = FALSE,
  sampling.f = c(.039, .027, .446, .008, 1e-6, .833, rep(1e-6, 10)),
  control = list(backend = "gslode")
)
st <- asr.marginal.classe(
  lik.td, mle.td$par.full,
  root = diversitree::ROOT.GIVEN, root.p = c(1, rep(0, 15))
)
rownames(st) <- state.labels

# Then extract node 495 and its matching full-tree ancestral reconstruction.
sisyphini <- extract.clade.asr(phy, node = 495, asr = st, tip.states = geo_regions)

# Keep five observed taxa, retaining full-tree ASR at surviving nodes.
sisyphini <- keep.tip.asr(sisyphini, tips = c(
  "Nesosisyphus_pygmaeus _STL3",
  "Nesosisyphus_regnardi _STL2",
  "Nesosisyphus_vicinus _STL1",
  "Nesosisyphus _rotundatus_STL52",
  "Sisyphus_sordidus_STL6"
))

sims <- make.simmap.classe.td(
  tree = sisyphini$tree, pars = mle.td$par.full,
  tip.weights = sisyphini$tip.states, node.probs = sisyphini$node.probs,
  k = 16, state.labels = state.labels, nsim = 100, seed = NULL
)
sim <- if (!is.null(sims$tree)) sims else sims[[1]]

# colors <- setNames(hcl.colors(length(state.labels), "Dynamic"), state.labels)
colors <- setNames(rep("grey80", length(state.labels)), state.labels)
colors[c("A", "S", "R", "A.S", "A.R", "S.R")] <- c(
  "#0072B2", "#E69F00", "#009E73",
  "#CC79A7", "#D55E00", "#56B4E9"
)
phytools::plotSimmap(sim$tree, colors = colors, split.vertical = TRUE, ftype = "i", fsize = 0.5)
#table(sim$tips$fate)

# Focus on the four Nesosisyphus paths and all side lineages born on them.
# Excludes the other root branch and every lineage born on that branch.
focus <- keep.lineage.simmap(sim, tips = c(
  "Nesosisyphus_pygmaeus _STL3",
  "Nesosisyphus_regnardi _STL2",
  "Nesosisyphus_vicinus _STL1",
  "Nesosisyphus _rotundatus_STL52"
))
phytools::plotSimmap(focus$tree, colors = colors, split.vertical = TRUE,
                    ftype = "i", fsize = 0.5)
legend(
  "topleft",
  legend = c("A", "S", "R", "A.S", "A.R", "S.R", "Other"),
  fill = c(colors[c("A", "S", "R", "A.S", "A.R", "S.R")], "grey80"),
  cex = 0.7, bty = "n"
)
#axisPhylo()
table(focus$tips$fate)

ltt <- lineages.through.time.simmap(
  sims,
  times = seq(9, 50, by = 1),
  states = c("S", "A.S", "S.R"),
  tips = focus$selected.tips
)

ltt$summary
plot(ltt)                          # Median + 95% simulation interval
plot(ltt, type = "maps")           # Translucent curve per map
plot(ltt, type = "segments",
     simulation = 1)              # Focal-state periods per branch
