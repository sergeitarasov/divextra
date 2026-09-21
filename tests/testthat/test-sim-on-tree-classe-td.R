test_that("backbone marking records ages, edges, and hard states", {
  b <- tiny_backbone()
  expect_s3_class(b, "classe_td_backbone")
  expect_equal(b$height, 0.5)
  expect_equal(b$age[b$root], 0.5)
  expect_equal(b$edge.table$parent.age - b$edge.table$child.age,
               b$edge.table$length)
  expect_equal(unname(b$node.states), 1L)
  expect_equal(unname(b$tip.states), c(1L, 1L))
  expect_false(b$probabilistic.nodes)
  expect_equal(unname(b$node.probs), matrix(1, 1, 1))
})

test_that("ancestral probability matrices remain soft node constraints", {
  tree <- ape::read.tree(text = "(a:1,(b:0.5,c:0.5):0.5);")
  nodes <- ape::Ntip(tree) + seq_len(tree$Nnode)
  node.probs <- matrix(
    c(1, 0, 0.2, 0.8), nrow = 2,
    dimnames = list(c("A", "B"), as.character(nodes))
  )
  backbone <- mark.classe.td.backbone(
    tree, node.probs, c(a = 1, b = 1, c = 1), 2, c("A", "B")
  )
  expect_true(backbone$probabilistic.nodes)
  expect_equal(backbone$node.probs[, as.character(nodes[2])], c(A = 0.2, B = 0.8))
  expect_equal(unname(backbone$node.states[as.character(nodes[2])]), 2L)

  pars <- c(
    lambda111.1 = 0.01, lambda112.1 = 0, lambda122.1 = 0,
    lambda211.1 = 0, lambda212.1 = 0, lambda222.1 = 0,
    mu1.1 = 0, mu2.1 = 0, q12.1 = 0, q21.1 = 0
  )
  sim <- simulate.classe.td.on.tree(
    backbone, pars, max.tries = 10, max.branch.tries = 100, seed = 4
  )
  event <- sim$events[
    sim$events$event == "observed_speciation" &
      sim$events$node == paste0("b", nodes[2]), , drop = FALSE
  ]
  expect_equal(event$from, 1L)
  expect_identical(sim$endpoint.constraint, "probability")
})

test_that("a supplied stem is mapped and restored after reconstruction", {
  tree <- ape::read.tree(text = "(a:0.5,b:0.5):0.4;")
  backbone <- mark.classe.td.backbone(
    tree, c(`3` = 1), c(a = 1, b = 1), 1, "A"
  )
  sim <- simulate.extinct.classe.td(
    backbone, tiny_pars(lambda = 0.2, mu = 2),
    max.tries = 10000, seed = 2
  )
  expect_equal(sum(sim$tree$root.map), sim$tree$root.edge)
  expect_equal(reconstruct.classe.td(sim)$tree$root.edge, tree$root.edge)
  expect_true(reconstruct.classe.td(sim)$matches.backbone)
})

test_that("epoch scheduling follows make.classe.td tip-to-root order", {
  schedule <- divextra:::.classe_td_schedule(tiny_pars(), 1, "A")
  expect_equal(divextra:::.classe_td_epoch(0.1, schedule$boundaries), 1L)
  expect_equal(divextra:::.classe_td_epoch(0.25, schedule$boundaries), 1L)
  expect_equal(divextra:::.classe_td_epoch(0.4, schedule$boundaries), 2L)
})

test_that("unconditioned candidate is a complete stochastic map", {
  sim <- simulate.classe.td.on.tree(
    tiny_backbone(), tiny_pars(lambda = 1.2, mu = 0.1),
    max.tries = 100, max.events = 10000, seed = 1
  )
  expect_valid_simmap(sim$tree, "A")
  expect_gt(sim$extra.extant, 0)
  expect_true(all(c("tip.status", "cladogenetic.events",
                    "event.history", "lineage.history") %in% names(sim$tree)))
  check <- reconstruct.classe.td(sim)
  expect_false(check$matches.backbone)
  expect_gt(length(check$extra.extant), 0)
})

test_that("observed backbone survival is conditioned rather than resimulated", {
  pars <- tiny_pars(lambda = 1e-8, mu = 1e6)
  sim <- simulate.classe.td.on.tree(
    tiny_backbone(), pars,
    max.tries = 1, max.branch.tries = 1, seed = 1
  )
  expect_true(sim$success)
  expect_equal(sim$attempts, 1L)
  expect_false(any(sim$events$event == "extinction" &
                     !is.na(sim$events$backbone.edge)))
})

test_that("nsim returns reproducible independent simulation objects", {
  sims <- simulate.classe.td.on.tree(
    tiny_backbone(), tiny_pars(lambda = 0.2, mu = 0),
    nsim = 3, max.tries = 100, seed = 12
  )
  expect_s3_class(sims, "classe_td_tree_simulations")
  expect_length(sims, 3)
  expect_true(all(vapply(sims, inherits, logical(1), "classe_td_tree_sim")))
  expect_equal(vapply(sims, `[[`, integer(1), "simulation"), 1:3)
  expect_length(reconstruct.classe.td(sims), 3)
})

test_that("conditioned simulation rejects complete histories with survivors", {
  sim <- simulate.extinct.classe.td(
    tiny_backbone(1), tiny_pars(lambda = 1.2, mu = 3, boundary = 0.5),
    max.tries = 10000, max.events = 10000, min.extinct = 1, seed = 1
  )
  expect_valid_simmap(sim$tree, "A")
  expect_true(sim$conditioned)
  expect_equal(sim$extra.extant, 0L)
  expect_true(reconstruct.classe.td(sim)$matches.backbone)
  expect_true(any(sim$tips$fate == "extinct"))
  expect_false(ape::is.ultrametric(sim$tree))
  expect_true(all(sim$tips$label[sim$tips$fate == "observed"] %in% c("a", "b")))
})

test_that("conditioned nsim produces multiple reconstructable maps", {
  sims <- simulate.extinct.classe.td(
    tiny_backbone(0.5), tiny_pars(lambda = 0.2, mu = 2),
    nsim = 2, max.tries = 10000, seed = 15
  )
  expect_s3_class(sims, "classe_td_tree_simulations")
  expect_true(all(vapply(
    reconstruct.classe.td(sims), `[[`, logical(1), "matches.backbone"
  )))
})

test_that("invalid rates and missing states fail early", {
  b <- tiny_backbone()
  bad <- tiny_pars()
  bad["mu1.1"] <- -1
  expect_error(simulate.classe.td.on.tree(b, bad), "non-negative")
  expect_error(simulate.classe.td.on.tree(b, tiny_pars(), nsim = 0), "nsim")

  tree <- ape::read.tree(text = "(a:1,b:1);")
  expect_error(
    mark.classe.td.backbone(tree, c(`3` = 1), c(a = 1), 1, "A"),
    "tip.states is missing"
  )
})

test_that("phytools can render an extinct-history simmap", {
  skip_if_not_installed("phytools")
  sim <- simulate.extinct.classe.td(
    tiny_backbone(1), tiny_pars(lambda = 1.2, mu = 3, boundary = 0.5),
    max.tries = 10000, seed = 1
  )
  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  expect_silent(phytools::plotSimmap(sim$tree, ftype = "off"))
})
