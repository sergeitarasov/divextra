test_that("Mauritian Neosisyphus fixture works with fitted 16-state model", {
  skip_on_cran()
  phy.path <- testthat::test_path("..", "..", "sim-on-tree", "phy-271.rds")
  fit.path <- testthat::test_path(
    "..", "..", "sim-on-tree", "AfrSunken16_gr-6_r-3.rds"
  )
  skip_if(!file.exists(phy.path) || !file.exists(fit.path),
          "local Neosisyphus example files are unavailable")

  phy <- readRDS(phy.path)
  fit <- readRDS(fit.path)
  n.tip <- ape::Ntip(phy)
  descendants <- function(node) {
    out <- integer()
    queue <- node
    while (length(queue)) {
      children <- phy$edge[phy$edge[, 1] %in% queue, 2]
      out <- c(out, children)
      queue <- children[children > n.tip]
    }
    out
  }
  tips <- descendants(513)
  tips <- tips[tips <= n.tip]
  tree <- ape::drop.tip(phy, phy$tip.label[-tips])
  states <- c(
    "A", "E", "M", "U", "S", "R", "A.E", "A.M", "A.U",
    "A.S", "A.R", "E.M", "E.U", "E.S", "E.R", "S.R"
  )
  node.states <- stats::setNames(
    rep(6L, tree$Nnode), ape::Ntip(tree) + seq_len(tree$Nnode)
  )
  tip.states <- stats::setNames(rep(6L, ape::Ntip(tree)), tree$tip.label)

  sim <- simulate.extinct.classe.td(
    tree, fit$par.full, node.states, tip.states,
    k = 16, state.labels = states, max.tries = 10000,
    min.extinct = 1, seed = 10
  )
  expect_valid_simmap(sim$tree, states)
  expect_true(reconstruct.classe.td(sim)$matches.backbone)
  expect_equal(sim$extra.extant, 0L)
  expect_gte(sum(sim$tips$fate == "extinct"), 1L)
})
