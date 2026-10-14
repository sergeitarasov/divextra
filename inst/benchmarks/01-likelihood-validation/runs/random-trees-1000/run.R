# Run from this directory: Rscript run.R ----
library(divextra)

# Use the archived source from this checkout, not a possibly older installation.
implementation <- new.env(parent = asNamespace("divextra"))
sys.source("provenance/classe_td.R", envir = implementation)
sys.source("provenance/sim-classe-td.R", envir = implementation)

set.seed(20261014)
n.tests <- 1000L
threshold <- 1e-5
control <- list(tol = 1e-10, eps = 1e-8)
dir.create("results", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)
inputs <- vector("list", n.tests)
results <- list()
row <- 0L
started <- proc.time()["elapsed"]

# Generate random numerical test cases, not model-conditioned datasets ----
for (i in seq_len(n.tests)) {
  n <- sample(100:500, 1)
  tree <- ape::rcoal(n)
  age <- runif(1, 1, 30)
  tree$edge.length <- tree$edge.length * age / max(ape::node.depth.edgelength(tree))
  states <- sample(rep(1:3, length.out = n))
  names(states) <- tree$tip.label
  geo.states <- states
  geo.states[states == 3] <- 0L # GeoSSE: AB=0, A=1, B=2; ClaSSE: A=1, B=2, AB=3.
  geosse.pars <- exp(runif(7, log(.005), log(.3)))
  names(geosse.pars) <- diversitree:::default.argnames.geosse()
  mapped.pars <- implementation$pars.geosse2classe(geosse.pars)
  full.pars <- exp(runif(27, log(.001), log(.1)))
  names(full.pars) <- diversitree:::default.argnames.classe(3)
  boundary <- age / 2
  inputs[[i]] <- list(tree = tree, states = states, geosse.pars = geosse.pars,
    mapped.pars = mapped.pars, full.pars = full.pars, boundary = boundary)
}
saveRDS(inputs, "inputs.rds")

# No optimization: evaluate the same trees and parameters ----
for (i in seq_len(n.tests)) {
  z <- inputs[[i]]
  geo.states <- z$states
  geo.states[z$states == 3] <- 0L
  episodic <- implementation$make.classe.td(z$tree, z$states,
    k = 3, n.epoch = 2, control = control)
  geosse <- diversitree::make.geosse(z$tree, geo.states, control = control)
  classe <- diversitree::make.classe(z$tree, z$states, k = 3, control = control)
  for (survival in c(FALSE, TRUE)) {
    values <- rep(NA_real_, 4)
    errors <- rep("", 4)
    # Each evaluation is retained, including errors and nonfinite results.
    for (j in 1:4) {
      ans <- tryCatch({
        if (j == 1) episodic(c(z$boundary, z$mapped.pars, z$mapped.pars),
          root = diversitree::ROOT.OBS, condition.surv = survival)
        else if (j == 2) geosse(z$geosse.pars,
          root = diversitree::ROOT.OBS, condition.surv = survival)
        else if (j == 3) episodic(c(z$boundary, z$full.pars, z$full.pars),
          root = diversitree::ROOT.OBS, condition.surv = survival)
        else classe(z$full.pars,
          root = diversitree::ROOT.OBS, condition.surv = survival)
      }, error = function(e) e)
      if (inherits(ans, "error")) errors[j] <- conditionMessage(ans)
      else values[j] <- as.numeric(ans)
    }
    for (comparison in 1:2) {
      a <- c(1, 3)[comparison]
      b <- c(2, 4)[comparison]
      row <- row + 1L
      difference <- values[a] - values[b]
      results[[row]] <- data.frame(id = i, n.tips = ape::Ntip(z$tree),
        age = 2 * z$boundary, survival = survival,
        comparison = c("GeoSSE", "ClaSSE")[comparison],
        episodic = values[a], reference = values[b], difference = difference,
        pass = is.finite(difference) && abs(difference) <= threshold,
        error = paste(errors[c(a,b)], collapse = ";"))
    }
  }
  if (i %% 50 == 0) {
    cat(i, "/", n.tests, "cases; elapsed", proc.time()["elapsed"] - started, "s\n")
    saveRDS(results, "results/checkpoint.rds")
  }
}
results <- do.call(rbind, results)
write.csv(results, "results/likelihood-checks.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "provenance/session-info.txt")
files <- c("inputs.rds", "run.R", "summarise.R", "provenance/classe_td.R",
  "provenance/sim-classe-td.R", "results/likelihood-checks.csv")
write.csv(data.frame(file = files, md5 = unname(tools::md5sum(files))),
  "provenance/checksums.csv", row.names = FALSE)
writeLines(as.character(proc.time()["elapsed"] - started), "results/elapsed-seconds.txt")
source("summarise.R")
