# Benchmark candidate and rejection simulation on a small reproducible fixture.
# Run from the package root with:
# Rscript inst/benchmarks/benchmark-sim-on-tree-classe-td.R

if (file.exists("R/sim-on-tree-classe-td.R")) {
  pkgload::load_all(".", export_all = FALSE)
} else {
  library(divextra)
}

tree <- ape::read.tree(text = "(a:0.5,b:0.5);")
backbone <- mark.classe.td.backbone(
  tree, c(`3` = 1), c(a = 1, b = 1), k = 1, state.labels = "A"
)
candidate.pars <- c(
  t.1 = 0.25,
  lambda111.1 = 1.2, mu1.1 = 0.1,
  lambda111.2 = 1.2, mu1.2 = 0.1
)
conditioned.pars <- c(
  t.1 = 0.25,
  lambda111.1 = 0.8, mu1.1 = 2,
  lambda111.2 = 0.8, mu1.2 = 2
)

run.benchmark <- function(n, fun) {
  elapsed <- numeric(n)
  attempts <- integer(n)
  events <- integer(n)
  sizes <- numeric(n)
  for (i in seq_len(n)) {
    timing <- system.time(ans <- fun(i))
    elapsed[i] <- unname(timing["elapsed"])
    attempts[i] <- ans$attempts
    events[i] <- ans$n.events
    sizes[i] <- as.numeric(utils::object.size(ans))
  }
  data.frame(
    n = n,
    simulations.per.second = n / sum(elapsed),
    median.seconds = stats::median(elapsed),
    p95.seconds = unname(stats::quantile(elapsed, 0.95)),
    mean.attempts = mean(attempts),
    p95.attempts = unname(stats::quantile(attempts, 0.95)),
    mean.events = mean(events),
    mean.object.bytes = mean(sizes)
  )
}

candidate <- run.benchmark(100, function(seed) {
  simulate.classe.td.on.tree(
    backbone, candidate.pars, max.tries = 1000,
    max.events = 100000, seed = seed
  )
})
conditioned <- run.benchmark(100, function(seed) {
  simulate.extinct.classe.td(
    backbone, conditioned.pars, max.tries = 10000,
    max.events = 100000, seed = seed
  )
})

print(rbind(candidate = candidate, conditioned = conditioned))
