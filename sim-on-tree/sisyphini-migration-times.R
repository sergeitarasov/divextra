# Uses the batch `sims` from sisyphini-stoch-mapping.R or
# sisyphini-lineages-through-time.R.
# Keep the four Nesosisyphus paths, then remove ALL generated side lineages.

focused <- lapply(sims, keep.lineage.simmap, tips = c(
  "Nesosisyphus_pygmaeus _STL3", "Nesosisyphus_regnardi _STL2",
  "Nesosisyphus_vicinus _STL1", "Nesosisyphus _rotundatus_STL52"
))
observed.maps <- lapply(focused, prune.unobserved.simmap)

phytools::plotSimmap(observed.maps[[3]]$tree, colors = colors,
                    split.vertical = TRUE, ftype = "i", fsize = .7)

migration <- transitions.through.time.simmap(
  observed.maps,
  breaks = seq(0, ceiling(max(ape::node.depth.edgelength(sisyphini$tree))), by = .5),
  from = "S.R",
  to = "R",
  event.types = c("anagenetic", "cladogenetic"), level = .95
)

migration$overall                # P(any arrival), MC CI and no-arrival map count
migration$age.summary            # First/last and pooled ages with 95% intervals
migration$event.age.distribution # Pooled-event age probabilities (event weighted)
migration$summary                # Per-bin counts and P(at least one arrival)
migration$per.map                # Includes separate counts for both event types
migration$events                 # Every arrival's age, from/to state and type

plot(migration, xlab = "Age (Ma)")
plot(migration, type = "first.age", xlab = "Age (Ma)")
plot(migration, type = "counts", xlab = "Age (Ma)")

# For A -> A.S alone, use from = "A", to = "A.S" above.
# For within-branch events alone, use event.types = "anagenetic".
