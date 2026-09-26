ltt_fixture <- function() {
  edge <- function(parent, child, start, end, map) list(parent = parent,
    child = child, start.age = start, end.age = end, length = start-end,
    map = map, backbone.edge = NA_integer_, generated = TRUE)
  edges <- list(edge("r", "a", 3, 0, c(S=1,A=1,S=1)),
                edge("r", "g", 3, 2, c(A=1)),
                edge("g", "b", 2, 0, c(S=2)),
                edge("g", "c", 2, 1, c(S=1)))
  tip <- function(node, age) data.frame(node=node, label=node, age=age,
    state=1L, fate=if(age==0) "extra_extant" else "extinct", generated=TRUE,
    state.label="S")
  .classe_td_graph_to_simmap(edges,
    list(a=tip("a",0),b=tip("b",0),c=tip("c",1)), c("S","A"))
}

test_that("counts use exact states, births, deaths and present-day tips", {
  sim <- ltt_fixture()
  times <- c(0,.5,1,1.5,2,2.5,3,4)
  x <- lineages.through.time.simmap(sim, times, "S")
  expect_equal(unname(x$counts[,1]), c(2,2,2,2,2,1,1,0))
  all <- lineages.through.time.simmap(sim, times, c("S","A"))
  expect_equal(unname(all$counts[,1]), c(2,2,2,3,3,2,2,0))
  expect_equal(x$summary$lower, x$summary$upper)
  expect_equal(sum(lineage.intervals.simmap(sim)$start.age -
                     lineage.intervals.simmap(sim)$end.age), sum(sim$tree$edge.length))
  expect_true(all(x$segments$state == "S"))
  expect_error(lineages.through.time.simmap(sim, times, "typo"), "Unknown")
  expect_error(lineages.through.time.simmap(sim, -1, "S"), "nonnegative")
  expect_error(lineages.through.time.simmap(sim, times, "S", probs=c(.9,.1)), "increasing")
})

test_that("zero maps remain in summaries and plots accept all three modes", {
  sim <- ltt_fixture()
  zero <- sim
  zero$tree$maps <- lapply(zero$tree$maps, function(m) setNames(sum(m), "A"))
  zero$tips$state.label <- "A"
  x <- lineages.through.time.simmap(list(sim,zero), c(0,1,2), "S")
  expect_equal(unname(x$counts[,2]), rep(0,3))
  expect_equal(x$summary$mean, c(1,1,1))
  expect_equal(x$summary$lower, rep(0,3))
  expect_equal(x$summary$upper, c(2,2,2))
  interpolated <- lineages.through.time.simmap(list(sim,zero),c(0,1,2),"S",quantile.type=7)
  expect_equal(interpolated$summary$lower,rep(.05,3))
  expect_equal(interpolated$summary$upper,rep(1.95,3))
  expect_equal(interpolated$counts,x$counts)
  expect_equal(interpolated$quantile.type,7)
  expect_error(lineages.through.time.simmap(sim,0:2,"S",quantile.type=1.5),"quantile.type")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(x))
  expect_invisible(plot(x, type="maps"))
  expect_invisible(plot(x, center="mean"))
  expect_invisible(plot(x, center="mean", level=.9))
  expect_invisible(plot(x, center="mean", interval="distribution"))
  expect_invisible(plot(x, center="mean", interval="distribution",quantile.type=7))
  expect_error(plot(x,quantile.type=NA),"quantile.type")
  expect_invisible(plot(x, type="maps", center="mean"))
  expect_error(plot(x, interval="mean-ci"), "requires center")
  expect_invisible(plot(x, type="segments"))
  expect_invisible(plot(x, type="segments", simulation=2))
})

test_that("batch focusing agrees with focusing each map in advance", {
  tree <- ape::read.tree(text="(a:1,b:1);")
  sims <- make.simmap.classe.td(tree, c(lambda111.1=.5,mu1.1=.2),
    tip.weights=c(a=1,b=1), node.probs=matrix(1,1,1,dimnames=list("1","3")),
    k=1, nsim=3, seed=2)
  a <- lineages.through.time.simmap(sims, seq(0,1,.1), "1", tips="a")
  b <- lineages.through.time.simmap(lapply(sims, keep.lineage.simmap, tips="a"),
                                   seq(0,1,.1), "1")
  expect_equal(a$counts, b$counts)
  expect_equal(a$segments, b$segments)
})

test_that("mean confidence intervals use Monte Carlo standard errors", {
  counts <- rbind(c(1,2,3,4),rep(2,4),rep(0,4))
  ci <- .classe_td_ltt_mean_ci(counts,.95)
  margin <- qt(.975,3)*sd(1:4)/2
  expect_equal(unname(ci[1,]),c(max(0,2.5-margin),2.5+margin))
  expect_equal(unname(ci[2,]),c(2,2))
  expect_equal(unname(ci[3,]),c(0,0))
  expect_lt(.classe_td_ltt_mean_ci(counts,.9)[1,2],ci[1,2])
  expect_error(.classe_td_ltt_mean_ci(counts,1),"level")
  expect_error(.classe_td_ltt_mean_ci(counts,NA_real_),"level")
  expect_error(.classe_td_ltt_mean_ci(counts[,1,drop=FALSE],.95),"at least two")
})
