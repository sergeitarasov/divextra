direct_extinction_schedule <- function(lambda = 0, mu = .7) {
  .classe_td_schedule(c(lambda111.1 = lambda, mu1.1 = mu), 1, "1")
}

test_that("shared extinction probabilities agree with analytic birth-death cases", {
  ages <- c(0, .001, .1, 1, 3)
  for (rates in list(c(0, .7), c(.3, .7), c(.7, .3), c(.5, .5))) {
    lambda <- rates[1]
    mu <- rates[2]
    cache <- .classe_td_extinction_cache(
      direct_extinction_schedule(lambda, mu), max.age = max(ages))
    expected <- if (lambda == mu) {
      lambda * ages / (1 + lambda * ages)
    } else {
      mu * (-expm1(-(lambda - mu) * ages)) /
        (lambda - mu * exp(-(lambda - mu) * ages))
    }
    actual <- vapply(ages, function(age) unname(cache$E(age)[1]), 0)
    expect_equal(actual, expected, tolerance = 2e-6)
    expect_true(all(actual >= 0 & actual <= 1))
  }
})

test_that("pure-death event ages follow a truncated exponential", {
  schedule <- direct_extinction_schedule(mu = .7)
  cache <- .classe_td_extinction_cache(schedule, max.age = 3)
  set.seed(824)
  waiting <- replicate(200, {
    event <- .classe_td_draw_extinct_event(3, 1, numeric(), schedule, cache)
    expect_identical(event$type, "extinction")
    expect_gt(event$age, 0)
    expect_lt(event$age, 3)
    expect_equal(unname(sum(event$map)), 3 - event$age, tolerance = 1e-8)
    3 - event$age
  })
  # Probability integral transform: a conservative fixed-seed ECDF bound.
  uniform <- -expm1(-.7 * waiting) / -expm1(-.7 * 3)
  expect_lt(max(abs(sort(uniform) - (seq_along(uniform) - .5) /
                        length(uniform))), .12)
})

test_that("piecewise extinction cache respects absolute epoch boundaries", {
  schedule <- .classe_td_schedule(
    c(lambda111.1 = 0, mu1.1 = .2,
      lambda111.2 = 0, mu1.2 = 1.1, t.1 = 1), 1, "1")
  cache <- .classe_td_extinction_cache(schedule, max.age = 3)
  ages <- c(0, .25, 1 - 1e-7, 1, 1 + 1e-7, 2, 3)
  expected <- -expm1(-(.2 * pmin(ages, 1) + 1.1 * pmax(ages - 1, 0)))
  expect_equal(vapply(ages, function(age) unname(cache$E(age)[1]), 0),
               expected, tolerance = 2e-6)
  set.seed(15)
  events <- replicate(100, .classe_td_draw_extinct_event(
    3, 1, numeric(), schedule, cache), simplify = FALSE)
  waits <- vapply(events, function(event) 3 - event$age, 0)
  # Conditional probability of reaching the epoch boundary without an event.
  expected.crossing <- exp(-1.1 * 2) * (-expm1(-.2)) / expected[length(expected)]
  expect_lt(abs(mean(waits > 2) - expected.crossing), .1)
  for (event in events) {
    expect_equal(event$epoch, .classe_td_epoch(event$age, schedule$boundaries))
    expect_equal(sum(event$map), 3 - event$age, tolerance = 1e-8)
  }
})

test_that("zero probability extinction errors but rare extinction is sampled", {
  schedule <- direct_extinction_schedule(lambda = .5, mu = 0)
  cache <- .classe_td_extinction_cache(schedule, max.age = 2)
  expect_equal(unname(cache$E(2)), 0)
  expect_error(.classe_td_draw_extinct_event(2, 1, numeric(), schedule, cache))
  rare <- direct_extinction_schedule(mu = 1e-8)
  rare.cache <- .classe_td_extinction_cache(rare, max.age = 2)
  expect_gt(rare.cache$E(2)[1], 0)
  set.seed(94)
  event <- .classe_td_draw_extinct_event(2, 1, numeric(), rare, rare.cache)
  expect_identical(event$type, "extinction")
  expect_true(event$age > 0 && event$age < 2)
})

test_that("state-dependent extinction agrees with an absorbing CTMC", {
  pars <- c(lambda111.1 = 0, lambda112.1 = 0, lambda122.1 = 0,
            lambda211.1 = 0, lambda212.1 = 0, lambda222.1 = 0,
            mu1.1 = 0, mu2.1 = 1.2, q12.1 = .8, q21.1 = .3)
  schedule <- .classe_td_schedule(pars, 2, c("A", "B"))
  cache <- .classe_td_extinction_cache(schedule, max.age = 3)
  transient <- matrix(c(-.8, .8, .3, -1.5), 2, 2, byrow = TRUE)
  for (age in c(.01, .1, 1, 3)) {
    expected <- 1 - rowSums(expm::expm(transient * age))
    expect_equal(unname(cache$E(age)), unname(expected), tolerance = 2e-6)
  }
  # A cannot die immediately, but conditioning is possible through A -> B.
  set.seed(51)
  event <- .classe_td_draw_extinct_event(3, 1, numeric(), schedule, cache)
  expect_identical(event$type, "transition")
  expect_equal(unname(event$to), 2)
})

direct_extinction_map_args <- function() {
  list(tree = ape::read.tree(text = "(a:2,b:2);"),
       pars = c(lambda111.1 = .5, mu1.1 = 1),
       tip.weights = c(a = 1, b = 1),
       node.probs = matrix(1, 1, 1, dimnames = list("1", "3")),
       k = 1)
}

test_that("E matches compiled ClaSSE likelihood across epochs and off-diagonal births", {
  nm <- diversitree:::default.argnames.classe(3)
  p1 <- setNames(rep(.03,length(nm)),nm)
  p1[c("mu1","mu2","mu3")] <- c(.2,.4,.6)
  p1[c("lambda112","lambda123","lambda213","lambda323")] <- c(.31,.19,.27,.43)
  p2 <- p1 * 1.7
  p2[c("mu1","mu2","mu3")] <- c(.5,.1,.3)
  pars <- c(t.1=2,setNames(p1,paste0(nm,".1")),setNames(p2,paste0(nm,".2")))
  tr <- ape::read.tree(text="((a:1,b:1):4,c:5);")
  lik <- make.classe.td(tr,c(a=1,b=2,c=3),k=3,n.epoch=2,
    sampling.f=rep(1,3),control=list(tol=1e-10))
  env <- environment(lik)
  branches <- environment(env$all_branches)$branches.td
  times <- c(.1,1,1.999,2,2.001,3,5)
  reference <- branches(c(rep(0,3),rep(1/3,3)),times,
    env$f.pars(pars),0,1L)[[2]][1:3,,drop=FALSE]
  schedule <- .classe_td_schedule(pars,3,c("A","B","AB"))
  cache <- .classe_td_extinction_cache(schedule,5)
  actual <- vapply(times,cache$E,numeric(3))
  expect_lt(max(abs(actual-reference)),1e-6)
  expect_gt(cache$E(2)[2],cache$E(5)[2]) # E itself need not be monotone.
  fine <- .classe_td_extinction_cache(schedule,5,list(max.step=.005,
    rtol=1e-11,atol=1e-14,interpolation.tol=1e-8))
  expect_lt(max(abs(vapply(times,fine$E,numeric(3))-reference)),
            max(abs(actual-reference)))
  # Conditioning at a fixed start must draw no survivors even across boundary2.
  for (seed in 1:20) {
    set.seed(seed)
    a <- .classe_td_draw_extinct_event(5,2,numeric(),schedule,cache)
    set.seed(seed)
    b <- .classe_td_draw_extinct_event(5,2,numeric(),schedule,fine)
    expect_identical(a$type,b$type)
    expect_equal(a$age,b$age,tolerance=1e-4)
  }
})

test_that("numerical controls cannot remove defaults or introduce invalid settings", {
  args <- c(direct_extinction_map_args(),list(side.lineage.mode="extinct_clades"))
  for (control in list(list(atol=NULL),list(max.step=-1),list(typo=1),
                      list(rtol=1e-8,rtol=1e-9)))
    expect_error(do.call(make.simmap.classe.td,c(args,list(extinction.control=control))),
                 "control|Control")
})

test_that("sampler selection leaves unrestricted mapping and RNG unchanged", {
  args <- c(direct_extinction_map_args(), list(seed = 32, nsim = 2))
  direct <- do.call(make.simmap.classe.td,
                   c(args, list(extinction.sampler = "direct")))
  direct.rng <- .Random.seed
  rejection <- do.call(make.simmap.classe.td,
                      c(args, list(extinction.sampler = "rejection")))
  expect_identical(.Random.seed, direct.rng)
  for (i in seq_along(direct)) {
    expect_identical(direct[[i]]$tree, rejection[[i]]$tree)
    expect_identical(direct[[i]]$events, rejection[[i]]$events)
    expect_identical(direct[[i]]$branch.attempts, rejection[[i]]$branch.attempts)
  }
})

test_that("direct and rejection side trees have compatible count distributions", {
  args <- c(direct_extinction_map_args(),
            list(side.lineage.mode = "extinct_clades", nsim = 100))
  direct <- do.call(make.simmap.classe.td,
    c(args, list(extinction.sampler = "direct", seed = 710)))
  rejection <- do.call(make.simmap.classe.td,
    c(args, list(extinction.sampler = "rejection", seed = 910)))
  counts <- function(sims) vapply(sims, function(sim) sum(sim$tips$generated), 0)
  a <- counts(direct)
  b <- counts(rejection)
  se <- sqrt(stats::var(a) / length(a) + stats::var(b) / length(b))
  expect_lt(abs(mean(a) - mean(b)), 5 * se + .05)
  for (sim in direct) {
    generated <- sim$tips$generated
    expect_true(all(sim$tips$fate[generated] == "extinct"))
    expect_true(all(sim$tips$age[generated] > 0))
    expect_equal(sim$extra.extant, 0)
    expect_equal(vapply(sim$tree$maps, sum, 0), sim$tree$edge.length)
    expect_setequal(sim$tips$label[!generated], c("a", "b"))
  }
})
