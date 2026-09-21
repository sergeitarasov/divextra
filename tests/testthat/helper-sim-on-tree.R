tiny_backbone <- function(length = 0.5) {
  tree <- ape::read.tree(text = sprintf("(a:%s,b:%s);", length, length))
  mark.classe.td.backbone(
    tree = tree,
    node.states = c(`3` = 1),
    tip.states = c(a = 1, b = 1),
    k = 1,
    state.labels = "A"
  )
}

tiny_pars <- function(lambda = 0.8, mu = 2, boundary = 0.25) {
  c(
    t.1 = boundary,
    lambda111.1 = lambda, mu1.1 = mu,
    lambda111.2 = lambda, mu1.2 = mu
  )
}

expect_valid_simmap <- function(tree, states) {
  expect_s3_class(tree, "classe_td_simmap")
  expect_true(inherits(tree, "simmap"))
  expect_true(inherits(tree, "phylo"))
  expect_length(tree$maps, nrow(tree$edge))
  expect_equal(length(tree$edge.length), nrow(tree$edge))
  expect_true(all(vapply(tree$maps, function(x) {
    is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
      all(x > 0) && !is.null(names(x)) && all(names(x) %in% states)
  }, logical(1))))
  expect_equal(vapply(tree$maps, sum, numeric(1)), tree$edge.length,
               tolerance = 1e-9)
  expect_equal(dim(tree$mapped.edge), c(nrow(tree$edge), length(states)))
  expect_equal(unname(rowSums(tree$mapped.edge)), tree$edge.length,
               tolerance = 1e-9)
}
