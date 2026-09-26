side_tip_args <- function() {
  tree <- ape::read.tree(text = "(((a:1,b:1):1,(c:1,d:1):1):1,e:3);")
  list(tree = tree, pars = c(lambda111.1 = .7, mu1.1 = .3),
       tip.weights = setNames(rep(1, 5), tree$tip.label),
       node.probs = matrix(1, 1, 4, dimnames = list("1", as.character(6:9))),
       k = 1, seed = 10)
}

test_that("all-tip selection preserves complete maps and RNG consumption", {
  args <- side_tip_args()
  args$nsim <- 2
  full <- do.call(make.simmap.classe.td, args)
  rng <- .Random.seed
  args$side.lineage.tips <- args$tree$tip.label
  all <- do.call(make.simmap.classe.td, args)
  expect_identical(.Random.seed, rng)
  for (i in seq_along(full)) {
    expect_identical(full[[i]]$tree, all[[i]]$tree)
    expect_identical(full[[i]]$branch.attempts, all[[i]]$branch.attempts)
    expect_true(all(all[[i]]$side.lineage.events$expanded))
  }
})

test_that("four-tip union retains shared stem once and excludes sister side trees", {
  args <- side_tip_args()
  full <- do.call(make.simmap.classe.td, args)
  args$side.lineage.tips <- c("a", "b", "c", "d")
  restricted <- do.call(make.simmap.classe.td, args)
  expect_equal(restricted$side.lineage.edges, which(args$tree$edge[, 2] != 5))
  expect_identical(restricted$branch.attempts, full$branch.attempts)
  # All observed nodes (including the root) and hidden backbone births remain.
  ev <- function(s) s$events[!is.na(s$events$backbone.edge) |
    s$events$event == "observed_speciation", setdiff(names(s$events), "phylo.node"), drop = FALSE]
  expect_equal(ev(restricted), ev(full))
  expect_setequal(restricted$tips$label[!restricted$tips$generated], args$tree$tip.label)
  births <- restricted$side.lineage.events
  expect_true(any(!births$expanded))
  roots <- restricted$lineages$parent[restricted$lineages$generated]
  expect_true(all(births$node[births$expanded] %in% roots))
  expect_false(any(births$node[!births$expanded] %in% roots))
  expect_equal(vapply(restricted$tree$maps, sum, 0), restricted$tree$edge.length)
  # Pruning all generated taxa recovers exactly the same mapped backbone.
  expect_equal(prune.unobserved.simmap(restricted)$tree,
               prune.unobserved.simmap(full)$tree)
  expect_error(do.call(make.simmap.classe.td, c(side_tip_args(),
    list(side.lineage.tips = "typo"))), "observed tip labels")
  expect_error(do.call(make.simmap.classe.td, c(side_tip_args(),
    list(side.lineage.tips = character()))), "observed tip labels")
})

test_that("restricted expansion still draws a state-changing root speciation", {
  p <- c(lambda111.1 = 0, lambda112.1 = 0, lambda122.1 = 1,
         lambda211.1 = 0, lambda212.1 = 0, lambda222.1 = .5,
         mu1.1 = 0, mu2.1 = .3, q12.1 = 0, q21.1 = 0)
  sim <- make.simmap.classe.td(ape::read.tree(text = "(a:1,b:1);"), p,
    tip.weights = c(a=2,b=2), node.probs = matrix(c(1,0),2,1,
    dimnames=list(c("A","R"),"3")), k=2, state.labels=c("A","R"),
    side.lineage.tips="a", seed=3)
  root <- sim$events[sim$events$event=="observed_speciation", ]
  expect_equal(root$from.label,"A")
  expect_equal(root$to1.label,"R")
  expect_equal(root$to2.label,"R")
  expect_equal(sim$tips$state[!sim$tips$generated],c(2L,2L))
})

test_that("focused counts agree distributionally with full expansion", {
  tree <- ape::read.tree(text = "(a:1,b:1);")
  args <- list(tree=tree, pars=c(lambda111.1=.5,mu1.1=.3),
    tip.weights=c(a=1,b=1), node.probs=matrix(1,1,1,dimnames=list("1","3")),k=1)
  # Select the SECOND daughter so skipped trees change the side-tree RNG draws.
  counts <- vapply(1:160, function(seed) {
    full <- do.call(make.simmap.classe.td,c(args,list(seed=seed)))
    part <- do.call(make.simmap.classe.td,c(args,list(seed=seed,side.lineage.tips="b")))
    f <- lineages.through.time.simmap(full,c(0,.5),"1",tips="b")$counts[,1]
    p <- lineages.through.time.simmap(part,c(0,.5),"1",tips="b")$counts[,1]
    f-p
  },numeric(2))
  se <- apply(counts,1,stats::sd)/sqrt(ncol(counts))
  expect_true(all(abs(rowMeans(counts)) < 5*se + .03))
})
