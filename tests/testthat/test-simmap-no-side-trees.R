test_that("skipping side trees preserves observed histories and conditioning", {
  tree <- ape::read.tree(text = "(a:2,b:2);")
  hidden <- 0L
  for (seed in 1:5) {
    args <- list(tree = tree, pars = c(lambda111.1 = .5, mu1.1 = 1),
      tip.weights = c(a = 1, b = 1), k = 1, state.labels = "A",
      backbone.method = "tree-conditioned", root.prior = 1, seed = seed)
    full <- do.call(make.simmap.classe.td, args)
    bare <- do.call(make.simmap.classe.td, c(args, list(expand.side.lineages = FALSE)))
    observed.full <- prune.unobserved.simmap(full)
    observed.bare <- prune.unobserved.simmap(bare)
    expect_equal(observed.bare$tree, observed.full$tree)
    expect_equal(bare$root.loglik, full$root.loglik)
    expect_equal(bare$root.p, full$root.p)
    expect_false(any(bare$tips$generated))
    expect_false(any(bare$side.lineage.events$expanded))
    hidden <- hidden + nrow(bare$side.lineage.events)
  }
  expect_gt(hidden, 0L)
})

test_that("no-side option preserves two-epoch GeoSSE backbone histories", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  states <- c(a = 1, b = 2, c = 3)
  old <- pars.geosse2classe(c(.12,.16,.08,.025,.035,.06,.09))
  recent <- pars.geosse2classe(c(.21,.18,.11,.045,.055,.13,.15))
  lik <- make.classe.td(tree, states, k = 3, n.epoch = 2)
  pars <- setNames(c(.7, recent, old), diversitree::argnames(lik))
  for (seed in 10:12) {
    args <- list(tree = tree, pars = pars, tip.weights = states, k = 3,
      state.labels = c("A","B","AB"), backbone.method = "tree-conditioned",
      root.prior = rep(1/3,3), seed = seed)
    full <- do.call(make.simmap.classe.td, args)
    bare <- do.call(make.simmap.classe.td, c(args, list(expand.side.lineages = FALSE)))
    expect_equal(prune.unobserved.simmap(full)$tree, prune.unobserved.simmap(bare)$tree)
    expect_equal(full$root.loglik, bare$root.loglik)
    expect_false(any(bare$tips$generated))
  }
  expect_error(do.call(make.simmap.classe.td, c(args, list(expand.side.lineages = NA))),
    "must be TRUE or FALSE")
})
