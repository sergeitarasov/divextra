# Representative likelihood comparisons before submission ----
source("model.R")
manifest <- read.csv("manifest.csv")
diag <- read.csv("tree-diagnostics.csv")
indices <- seq(1, 750, by = 50)
indices <- unique(c(indices, diag$tree.index[diag$boundary.before.crown]))
checks <- data.frame()
for (i in indices) {
  z <- readRDS(sprintf("datasets/tree-%03d.rds", i))
  stopifnot(ape::Ntip(z$tree) == manifest$N[i])
  p <- z$truth[1:3]
  geosse <- make.likelihood(z, "geosse")
  reference <- geosse(unname(p))
  for (fraction in c(.05, .25, .5)) {
    z$recent.time <- z$duration * fraction
    linked <- make.likelihood(z, "linked_classe")
    ll <- linked(p)
    checks <- rbind(checks, data.frame(tree.index = i, fraction = fraction,
      geosse = reference, classe = ll, difference = ll - reference))
  }
}
write.csv(checks, "preparation/likelihood-checks.csv", row.names = FALSE)
stopifnot(all(is.finite(checks$classe)), max(abs(checks$difference)) < 1e-4)
cat("Likelihoods checked on", length(indices), "trees; maximum difference:", max(abs(checks$difference)), "\n")

# One actual optimization using the same fit function as Roihu ----
z <- readRDS("datasets/tree-001.rds")
job <- list(model = "geosse", start = 1L, initial = c(s = .08, d = .08, x = .03))
fit <- fit.one(z, job)
saveRDS(fit, "preparation/local-test-fit.rds")
stopifnot(is.finite(fit$lnLik), fit$lnLik > -1e49)
print(fit[c("par", "lnLik", "convergence", "elapsed")])
