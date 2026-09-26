geosse_td_example_parameters <- function() {
  pars.ge <- matrix(
    c(
      0.1, 0.1, 0.1, 0, 0, 0.1, 0.1,
      0.3, 0.3, 0.3, 0, 0, 0.1, 0.1
    ),
    2, 7, byrow = TRUE
  )
  colnames(pars.ge) <- diversitree:::default.argnames.geosse()
  t(apply(pars.ge, 1, pars.geosse2classe))
}

test_that("table2simmap preserves simulated transition histories", {
  set.seed(123)
  info <- make.tree.classe.td(
    geosse_td_example_parameters(), k = 3,
    max.t1 = 10, max.t2 = 15, x0 = 1, single.lineage = TRUE
  )
  expect_gt(nrow(attr(info, "hist")), 0)

  labels <- c("A", "B", "AB")
  simmap <- table2simmap(info, state.labels = labels)

  expect_s3_class(simmap, "simmap")
  expect_s3_class(simmap, "phylo")
  expect_length(simmap$maps, nrow(simmap$edge))
  expect_equal(
    unname(vapply(simmap$maps, sum, numeric(1))),
    unname(simmap$edge.length),
    tolerance = 1e-10
  )
  expect_equal(
    unname(rowSums(simmap$mapped.edge)),
    unname(simmap$edge.length),
    tolerance = 1e-10
  )
  expect_identical(colnames(simmap$mapped.edge), labels)
  expect_gt(sum(vapply(simmap$maps, length, integer(1)) > 1L), 0)

  tip.edges <- which(simmap$edge[, 2] <= ape::Ntip(simmap))
  terminal.map.states <- vapply(
    simmap$maps[tip.edges], function(map) names(map)[length(map)], character(1)
  )
  expect_identical(
    unname(terminal.map.states),
    unname(labels[
      simmap$tip.state[simmap$tip.label][simmap$edge[tip.edges, 2]]
    ])
  )
})

test_that("table2simmap can be displayed by phytools", {
  skip_if_not_installed("phytools")
  set.seed(123)
  info <- make.tree.classe.td(
    geosse_td_example_parameters(), k = 3,
    max.t1 = 10, max.t2 = 15, x0 = 1, single.lineage = TRUE
  )
  simmap <- table2simmap(info, state.labels = c("A", "B", "AB"))
  colors <- c(A = "#D55E00", B = "#0072B2", AB = "#009E73")

  path <- tempfile(fileext = ".pdf")
  grDevices::pdf(path)
  on.exit({
    grDevices::dev.off()
    unlink(path)
  }, add = TRUE)
  expect_silent(phytools::plotSimmap(simmap, colors = colors, ftype = "off"))
})

test_that("table2simmap retains maps when extinct tips are pruned", {
  pars.ge <- matrix(
    c(
      0.15, 0.15, 0.15, 0.08, 0.08, 0.1, 0.1,
      0.30, 0.30, 0.30, 0.08, 0.08, 0.1, 0.1
    ),
    2, 7, byrow = TRUE
  )
  colnames(pars.ge) <- diversitree:::default.argnames.geosse()
  pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))

  set.seed(6)
  info <- make.tree.classe.td(
    pars.cl, k = 3, max.t1 = 5, max.t2 = 10,
    x0 = 1, single.lineage = TRUE
  )
  expect_gt(sum(!info$split & info$extinct), 0)

  simmap <- table2simmap(info, state.labels = c("A", "B", "AB"))
  expect_false(any(startsWith(simmap$tip.label, "ex")))
  expect_equal(
    unname(vapply(simmap$maps, sum, numeric(1))),
    unname(simmap$edge.length),
    tolerance = 1e-10
  )
  # One retained edge crosses both a simulated transition and a suppressed
  # cladogenetic event, so all three chronological segments must survive.
  expect_true(any(vapply(simmap$maps, length, integer(1)) >= 3L))
})
