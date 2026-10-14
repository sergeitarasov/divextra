# Diagnostic rerun only: never replaces the primary benchmark results ----
library(divextra)
implementation <- new.env(parent = asNamespace("divextra"))
sys.source("provenance/classe_td.R", envir = implementation)
inputs <- readRDS("inputs.rds")
checks <- read.csv("results/likelihood-checks.csv")
failed <- unique(checks$id[!checks$pass])
rechecks <- data.frame()
for (i in failed) {
  z <- inputs[[i]]
  states <- z$states
  states[states == 3] <- 0L
  episodic <- implementation$make.classe.td(z$tree, z$states, k = 3,
    n.epoch = 2, control = list(tol = 1e-10, eps = 1e-8))
  geosse <- diversitree::make.geosse(z$tree, states,
    control = list(tol = 1e-10, eps = 1e-8))
  for (survival in c(FALSE, TRUE)) {
    a <- tryCatch(episodic(c(z$boundary, z$mapped.pars, z$mapped.pars),
      condition.surv = survival, root = diversitree::ROOT.OBS), error = conditionMessage)
    b <- tryCatch(geosse(z$geosse.pars, condition.surv = survival,
      root = diversitree::ROOT.OBS), error = conditionMessage)
    rechecks <- rbind(rechecks, data.frame(id = i, survival = survival,
      episodic = as.character(a), reference = as.character(b)))
  }
}
write.csv(rechecks, "results/error-rechecks.csv", row.names = FALSE)
print(rechecks)
