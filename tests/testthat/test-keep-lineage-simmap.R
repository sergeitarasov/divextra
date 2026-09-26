test_that("lineage filtering excludes other root branch and all its births", {
  tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
  sim <- make.simmap.classe.td(tree,
    pars = c(lambda111.1 = 1, mu1.1 = .5),
    tip.weights = setNames(rep(1, 3), tree$tip.label),
    node.probs = matrix(1, 1, 2, dimnames = list("1", c("4", "5"))),
    k = 1, seed = 17)
  out <- keep.lineage.simmap(sim, c("a", "b"))
  lin <- sim$lineages
  # Independently collect everything below the unwanted backbone root child.
  edge.c <- which(tree$edge[, 2] == 3)
  bad.nodes <- lin$child[!lin$generated & lin$backbone.edge == edge.c &
                         !is.na(lin$backbone.edge)]
  repeat {
    next.nodes <- unique(c(bad.nodes, lin$child[lin$parent %in% bad.nodes]))
    if (setequal(next.nodes, bad.nodes)) break
    bad.nodes <- next.nodes
  }
  expect_true(any(sim$tips$generated & sim$tips$node %in% bad.nodes))
  expect_setequal(out$tips$label, sim$tips$label[!sim$tips$node %in% bad.nodes])
  expect_equal(out$tips$label[!out$tips$generated], c("a", "b"))
  expect_equal(vapply(out$tree$maps, sum, 0), out$tree$edge.length)
  expect_true(all(out$lineages$start.age %in% sim$lineages$start.age))
  expect_equal(out$lineages$map,
    lin$map[match(out$lineages$child, lin$child)])
  expect_setequal(unname(keep.lineage.simmap(sim, tree$tip.label)$tree$tip.label),
                  unname(sim$tree$tip.label))
  one <- keep.lineage.simmap(sim, "a")
  expect_equal(one$tips$label[!one$tips$generated], "a")
  expect_error(keep.lineage.simmap(sim, "unknown"), "observed backbone")
})

test_that("batches match individual focusing and preserve names and shape", {
  tree <- ape::read.tree(text = "(a:1,b:1);")
  sims <- make.simmap.classe.td(tree,
    pars = c(lambda111.1 = .5, mu1.1 = .2),
    tip.weights = c(a = 1, b = 1),
    node.probs = matrix(1, 1, 1, dimnames = list("1", "3")),
    k = 1, nsim = 2, seed = 10)
  names(sims) <- c("first", "second")
  focused <- keep.lineage.simmap(sims, "a")
  expect_identical(focused, lapply(sims, keep.lineage.simmap, tips = "a"))
  expect_identical(names(focused), names(sims))
  expect_length(keep.lineage.simmap(list(sims[[1]]), "a"), 1)
  expect_equal(keep.lineage.simmap(list(sims[[1]]), "a")[[1]], focused[[1]])
  expect_error(keep.lineage.simmap(list(), "a"), "nonempty")
  expect_error(keep.lineage.simmap(list(sims[[1]], NULL), "a"), "index 2")
  expect_error(keep.lineage.simmap(sims, "unknown"), "Simulation 1")
  expect_equal(lineages.through.time.simmap(focused, c(0,.5,1), "1")$counts,
    lineages.through.time.simmap(sims, c(0,.5,1), "1", tips = "a")$counts)
})

test_that("a single selected tip without side births retains its whole history", {
  tree <- ape::read.tree(text = "(a:1,b:1);")
  sim <- make.simmap.classe.td(tree,
    pars = c(lambda111.1 = 1e-12, mu1.1 = 0),
    tip.weights = c(a = 1, b = 1),
    node.probs = matrix(1, 1, 1, dimnames = list("1", "3")), k = 1, seed = 1)
  out <- keep.lineage.simmap(sim, "a")
  expect_equal(unname(out$tree$tip.label), "a")
  expect_equal(sum(out$tree$edge.length), 1)
  expect_equal(sum(out$tree$maps[[1]]), 1)
})
