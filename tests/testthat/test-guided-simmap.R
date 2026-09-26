guided_test_input <- function(k = 1L, length = 1) {
  tree <- ape::read.tree(text = sprintf("(a:%s,b:%s);", length, length))
  list(
    tree = tree,
    tip.weights = matrix(1, k, 2, dimnames = list(NULL, tree$tip.label)),
    node.probs = matrix(1 / k, k, 1, dimnames = list(NULL, "3")),
    root.p = c(1, rep(0, k - 1L)),
    k = k
  )
}

guided_test_two_state_pars <- function() {
  c(lambda111.1 = 0, lambda112.1 = 0, lambda122.1 = 0,
    lambda211.1 = 0, lambda212.1 = 0, lambda222.1 = 0,
    mu1.1 = 0, mu2.1 = 0, q12.1 = 0, q21.1 = 0)
}

test_that("guided mapping rejects an impossible observed root event", {
  args <- guided_test_input()
  args$pars <- c(lambda111.1 = 0, mu1.1 = 0)
  expect_error(do.call(make.simmap.classe.td, args),
               "[Ss]peciation|lambda|[Rr]oot|incompatible_observed_node")
})

test_that("explicit root probabilities override the root ASR column", {
  args <- guided_test_input(2)
  args$pars <- guided_test_two_state_pars()
  args$pars[c("lambda111.1", "lambda222.1")] <- 0.1
  args$node.probs[, 1] <- c(1, 0)
  args$root.p <- c(0, 1)
  sim <- do.call(make.simmap.classe.td, c(args, list(seed = 12)))
  root.event <- sim$events[sim$events$event == "observed_speciation", ]
  expect_equal(root.event$from, 2L)
  expect_true(all(vapply(sim$tree$maps, function(x) {
    all(names(x) == "2")
  }, logical(1))))
})

test_that("root ASR is used by default and custom root is reproducible", {
  args <- guided_test_input(2)
  args$root.p <- NULL
  args$node.probs[, 1] <- c(0, 1)
  args$pars <- guided_test_two_state_pars()
  args$pars[c("lambda111.1", "lambda222.1")] <- .1
  implicit <- do.call(make.simmap.classe.td, c(args, list(seed = 10)))
  expect_equal(implicit$root.state, 2L)
  expect_identical(implicit$root.source, "node-ASR")
  args$root.p <- c(0, 1)
  explicit <- do.call(make.simmap.classe.td, c(args, list(seed = 10)))
  expect_identical(explicit$root.source, "custom")
  expect_equal(implicit$tree, explicit$tree)
})

test_that("endpoint rejection redraws the entire daughter pair", {
  args <- guided_test_input(2)
  args$pars <- guided_test_two_state_pars()
  args$pars[c("lambda111.1", "lambda122.1")] <- 1e-12
  args$tip.weights[,] <- c(.2, .8)
  sims <- do.call(make.simmap.classe.td, c(args, list(nsim = 200, seed = 982)))
  pairs <- t(vapply(sims, function(s) {
    e <- s$events[s$events$event == "observed_speciation", ]
    as.integer(c(e$to1, e$to2))
  }, integer(2)))
  expect_true(all(pairs[, 1] == pairs[, 2]))
  expected <- .8^2 / (.2^2 + .8^2)
  expect_lt(abs(mean(pairs[, 1] == 2) - expected),
            5 * sqrt(expected * (1 - expected) / nrow(pairs)))
  expect_true(any(vapply(sims, function(s) sum(s$branch.attempts) > 1, logical(1))))
})

test_that("soft tip weights are accepted and sampled as endpoint weights", {
  args <- guided_test_input(2)
  args$pars <- guided_test_two_state_pars()
  args$pars["lambda111.1"] <- 1e-12
  args$pars[c("q12.1", "q21.1")] <- c(0.7, 0.2)
  args$tip.weights[,] <- c(0.2, 0.8)
  sims <- do.call(make.simmap.classe.td,
                 c(args, list(nsim = 250, seed = 318)))
  endpoint.two <- unlist(lapply(sims, function(sim) {
    terminal.edges <- match(match(c("a", "b"), sim$tree$tip.label),
                            sim$tree$edge[, 2])
    vapply(sim$tree$maps[terminal.edges], function(x) {
      tail(names(x), 1) == "2"
    }, logical(1))
  }))
  p12 <- 0.7 / 0.9 * (1 - exp(-0.9))
  expected <- p12 * 0.8 / ((1 - p12) * 0.2 + p12 * 0.8)
  # The target is K_ij*w_j normalized, not the supplied marginal 0.8.
  se <- sqrt(expected * (1 - expected) / length(endpoint.two))
  expect_lt(abs(mean(endpoint.two) - expected), 5 * se)
})

test_that("zero extinction preserves all generated lineages at present", {
  args <- guided_test_input()
  args$pars <- c(lambda111.1 = 1.2, mu1.1 = 0)
  sims <- do.call(make.simmap.classe.td,
                 c(args, list(nsim = 12, seed = 82)))
  expect_length(sims, 12)
  expect_true(all(vapply(sims, function(sim) {
    !any(sim$tips$fate == "extinct")
  }, logical(1))))
  expect_gt(sum(vapply(sims, function(sim) {
    sum(sim$tips$fate == "extra_extant")
  }, integer(1))), 0)
})

test_that("distinct daughters are assigned to either backbone side equally", {
  args <- guided_test_input(2)
  args$pars <- guided_test_two_state_pars()
  args$pars[c("lambda112.1", "lambda222.1")] <- 1e-12
  sims <- do.call(make.simmap.classe.td,
                 c(args, list(nsim = 160, seed = 441)))
  left.is.one <- vapply(sims, function(sim) {
    edge <- which(sim$tree$edge[, 2] == match("a", sim$tree$tip.label))
    names(sim$tree$maps[[edge]])[1] == "1"
  }, logical(1))
  expect_lt(abs(mean(left.is.one) - 0.5), 5 * sqrt(0.25 / length(sims)))
})

test_that("hidden speciation can replace the backbone parental state", {
  args <- guided_test_input(2)
  args$pars <- guided_test_two_state_pars()
  args$pars[c("lambda122.1", "lambda211.1")] <- 3
  args$pars[c("mu1.1", "mu2.1")] <- 6
  sims <- do.call(make.simmap.classe.td,
                 c(args, list(nsim = 5, seed = 173)))
  hidden <- do.call(rbind, lapply(sims, function(sim) {
    sim$events[sim$events$event == "hidden_speciation", , drop = FALSE]
  }))
  expect_gt(nrow(hidden), 0)
  # In this model every allowed split has two daughters unlike the parent.
  expect_true(all(hidden$from != hidden$to1 & hidden$from != hidden$to2))
  for (sim in sims) expect_valid_simmap(sim$tree, c("1", "2"))
})

test_that("unrestricted side trees include extinct and extra extant tips", {
  args <- guided_test_input()
  args$pars <- c(lambda111.1 = 1.2, mu1.1 = 1.5)
  sims <- do.call(make.simmap.classe.td,
                 c(args, list(nsim = 30, seed = 159)))
  fates <- unlist(lapply(sims, function(sim) sim$tips$fate))
  expect_true(all(c("extinct", "extra_extant", "observed") %in% fates))
  expect_true(all(vapply(sims, function(sim) {
    observed <- sim$tips$label[sim$tips$fate == "observed"]
    setequal(observed, args$tree$tip.label)
  }, logical(1))))
})

test_that("one-state backbone hidden births follow their integrated rate", {
  args <- guided_test_input(length = 0.5)
  args$pars <- c(t.1 = 0.25,
                 lambda111.1 = 1.6, mu1.1 = 4,
                 lambda111.2 = 0.8, mu1.2 = 4)
  sims <- do.call(make.simmap.classe.td,
                 c(args, list(nsim = 180, seed = 19)))
  counts <- vapply(sims, function(sim) {
    sum(sim$events$event == "hidden_speciation")
  }, integer(1))
  expected <- 2 * (0.25 * 1.6 + 0.25 * 0.8)
  expect_lt(abs(mean(counts) - expected), 5 * sqrt(expected / length(sims)))
})

test_that("complete maps preserve the supplied three-tip backbone", {
  args <- guided_test_input()
  args$tree <- ape::read.tree(text = "(a:1,(b:0.4,c:0.4):0.6);")
  args$tip.weights <- matrix(1, 1, 3,
                             dimnames = list(NULL, args$tree$tip.label))
  args$node.probs <- matrix(1, 1, 2,
                            dimnames = list(NULL, c("4", "5")))
  args$pars <- c(lambda111.1 = 1.2, mu1.1 = 1.5)
  sims <- do.call(make.simmap.classe.td,
                 c(args, list(nsim = 8, seed = 62)))
  for (sim in sims) {
    expect_valid_simmap(sim$tree, "1")
    # Drop all generated tips, including survivors, to recover the input tree.
    recovered <- ape::keep.tip(sim$tree, args$tree$tip.label)
    expect_equal(ape::cophenetic.phylo(recovered)[args$tree$tip.label,
                                                args$tree$tip.label],
                 ape::cophenetic.phylo(args$tree))
    expect_equal(unname(ape::dist.nodes(recovered)[ape::Ntip(recovered) + 1L,
                                          seq_len(ape::Ntip(recovered))]),
                 rep(1, 3))
  }
})

test_that("an impossible endpoint fails at the local retry limit", {
  args <- guided_test_input(2)
  args$pars <- guided_test_two_state_pars()
  args$pars["lambda111.1"] <- 0.1
  args$tip.weights[,] <- c(0, 1)
  expect_error(do.call(make.simmap.classe.td,
                       c(args, list(max.branch.tries = 2, seed = 77))),
               "[Rr]etr|[Aa]ttempt|[Ee]ndpoint|[Bb]ranch|tries")
})
