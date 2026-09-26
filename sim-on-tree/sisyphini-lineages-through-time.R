# First run sisyphini-stoch-mapping.R to prepare sisyphini, mle.td and state.labels.
# No files are saved. Use a fixed seed to reproduce this batch.
sims <- make.simmap.classe.td(
  tree = sisyphini$tree, pars = mle.td$par.full,
  tip.weights = sisyphini$tip.states, node.probs = sisyphini$node.probs,
  k = 16, state.labels = state.labels, nsim = 100, seed = 10
)

ltt <- lineages.through.time.simmap(
  sims,
  times = seq(0, max(ape::node.depth.edgelength(sisyphini$tree)), by = 1),
  states = c("S", "A.S", "S.R"),
  tips = c("Nesosisyphus_pygmaeus _STL3", "Nesosisyphus_regnardi _STL2",
           "Nesosisyphus_vicinus _STL1", "Nesosisyphus _rotundatus_STL52"),
  probs = c(.025, .975)
)

# Pooled S / A.S / S.R counts: median and pointwise 95% simulation interval.
ltt$summary
plot(ltt, xlab = "Age (Ma)")

# Each translucent curve represents the counts from one simulated map.
plot(ltt, type = "maps", xlab = "Age (Ma)")

# Exact periods spent in focal states on each branch lineage of map 1.
# A row is a branch segment between nodes, not species identity across nodes.
plot(ltt, type = "segments", simulation = 1, xlab = "Age (Ma)")
