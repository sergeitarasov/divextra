extinction_args <- function() {
  list(tree = ape::read.tree(text = "(a:6,b:6);"),
       pars = c(lambda111.1 = .4, mu1.1 = 1.5,
                lambda111.2 = .7, mu1.2 = 2, t.1 = 4),
       tip.weights = c(a=1,b=1),
       node.probs = matrix(1,1,1,dimnames=list("1","3")),
       k = 1, seed = 10)
}

test_that("default all mode preserves RNG and mapped histories", {
  args <- extinction_args()
  a <- do.call(make.simmap.classe.td, args)
  rng <- .Random.seed
  b <- do.call(make.simmap.classe.td, c(args,list(side.lineage.mode="all")))
  expect_identical(a, b)
  expect_identical(.Random.seed, rng)
})

test_that("side clades die before present across epochs and batches", {
  for (sampler in c("direct","rejection")) {
    sims <- do.call(make.simmap.classe.td, c(extinction_args(),
      list(side.lineage.mode="extinct_clades", extinction.sampler=sampler, nsim=8)))
    n.generated <- 0L
    for (sim in sims) {
      generated <- sim$tips[sim$tips$generated, , drop=FALSE]
      n.generated <- n.generated + nrow(generated)
      expect_true(all(generated$fate == "extinct"))
      expect_true(all(generated$age > 0))
      expect_true(all(sim$side.lineage.events$age > 0))
      expect_equal(sim$extra.extant, 0L)
      expect_equal(sim$extinction.sampler, sampler)
      expect_setequal(sim$tips$label[!sim$tips$generated], c("a","b"))
      expect_equal(sum(sim$events$event == "observed_speciation"), 1L)
      expect_equal(vapply(sim$tree$maps,sum,0), sim$tree$edge.length)
      expect_equal(length(sim$side.attempts), nrow(sim$side.lineage.events))
      expect_false(anyDuplicated(sim$tips$label) > 0L)
      expect_false(anyDuplicated(sim$events$event.id) > 0L)
    }
    expect_gt(n.generated, 0)
  }
})

test_that("zero deadline preserves backbone; selected paths expand only their sides", {
  args <- extinction_args()
  full <- do.call(make.simmap.classe.td, args)
  conditioned <- do.call(make.simmap.classe.td, c(args,
    list(side.lineage.mode="extinct_clades", side.lineage.tips="a")))
  expect_identical(full$branch.attempts, conditioned$branch.attempts)
  expect_equal(prune.unobserved.simmap(full)$tree,
               prune.unobserved.simmap(conditioned)$tree)
  expect_equal(length(conditioned$side.attempts),
               sum(conditioned$side.lineage.events$expanded))
})

test_that("impossible extinction and resource exhaustion fail explicitly", {
  args <- extinction_args()
  args$pars <- c(lambda111.1=2, mu1.1=0)
  expect_error(do.call(make.simmap.classe.td,c(args,
    list(side.lineage.mode="extinct_clades", extinction.sampler="rejection",
         max.side.tries=2))), "max_side_tries")
  expect_error(do.call(make.simmap.classe.td,c(args,
    list(side.lineage.mode="extinct_clades", max.events=1))), "max_events")
})

test_that("deadline removed and retry arguments are validated", {
  expect_error(do.call(make.simmap.classe.td,c(extinction_args(),
    list(extinction.age=3))), "unused argument")
  expect_error(do.call(make.simmap.classe.td,c(extinction_args(),
    list(max.side.tries=0))), "max.side.tries")
})
