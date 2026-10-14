# Read-only remote check: no optimization on the login node ----
source("model.R")
stopifnot(requireNamespace("callr",quietly=TRUE))
definitions <- readRDS("provenance/runtime-functions.rds")
for(name in names(definitions)) {
  f <- get(name,asNamespace("divextra"))
  stopifnot(identical(definitions[[name]],c(deparse(formals(f)),deparse(body(f)))))
}
checks <- read.csv("preparation/likelihood-checks.csv")
for(i in c(1,501,1001)) {
  z <- readRDS(sprintf("datasets/tree-%03d.rds",i))
  evaluate <- make.likelihood(z,"linked_classe")
  ll <- evaluate(z$truth[c("s_recent","d_recent","x_recent")])
  reference <- checks$lnLik[checks$tree.index==i & checks$perturbation==0]
  stopifnot(length(reference)==1,abs(ll-reference)<1e-5)
  cat(z$scenario,"remote minus local likelihood:",ll-reference,"\n")
}
