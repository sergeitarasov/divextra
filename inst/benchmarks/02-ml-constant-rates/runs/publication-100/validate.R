# Run from the benchmark directory after simulate.R ----
source("model.R")
m <- read.csv("manifest.csv")
stopifnot(nrow(m)==1500, !anyDuplicated(m$seed), all(table(m$scenario,m$tip.group)==100),
  all(m$N>=50 & m$N<=2000))
checks <- data.frame()
for(i in seq(1,1500,by=100)) {
  z <- readRDS(sprintf("datasets/tree-%03d.rds",i))
  p <- z$truth[c("s_recent","d_recent","x_recent")]
  evaluate <- make.likelihood(z,"linked_classe")
  # Linked constant rates must be invariant to a bookkeeping boundary shift.
  shifted <- z
  shifted$recent.time <- z$duration*.25
  evaluate.shifted <- make.likelihood(shifted,"linked_classe")
  truth.ll <- evaluate(p)
  for(j in 0:3) {
    altered <- p
    if(j>0) altered[j] <- altered[j]*1.1
    ll <- evaluate(altered)
    other <- evaluate.shifted(altered)
    stopifnot(is.finite(ll), abs(ll-other)<1e-5)
    if(j>0) stopifnot(abs(ll-truth.ll)>1e-8)
    checks <- rbind(checks,data.frame(tree.index=i,scenario=z$scenario,
      perturbation=j,lnLik=ll,shifted.lnLik=other))
  }
}
write.csv(checks,"preparation/likelihood-checks.csv",row.names=FALSE)

# Three starts on a small tree at each rate scale; same timeout as cluster ----
fits <- list()
for(i in c(1,501,1001)) {
  z <- readRDS(sprintf("datasets/tree-%03d.rds",i))
  for(start in 1:3) {
    set.seed(z$seed+start*10000000)
    initial <- c(s_recent=.08,d_recent=.08,x_recent=.03)*c(.5,1,2)[start]*
      exp(runif(3,-log(2),log(2)))
    job <- list(model="linked_classe",start=start,initial=initial)
    f <- callr::r(function(z,job) { source("model.R"); fit.one(z,job) },
      args=list(z=z,job=job),libpath=.libPaths(),timeout=1200)
    stopifnot(is.finite(f$lnLik),f$lnLik > -1e49,length(f$par)==3)
    fits[[paste(z$scenario,start)]] <- f
    print(f[c("id","par","lnLik","convergence")])
  }
}
saveRDS(fits,"preparation/local-three-starts.rds")

# Capture exact runtime function definitions, without private environment data.
functions <- c("make.tree.classe.td","make.classe.td","pars.geosse2classe","make_constraints_sse_td","read_yaml_pars_td")
definitions <- list()
for(name in functions) {
  f <- get(name,asNamespace("divextra"))
  definitions[[name]] <- c(deparse(formals(f)),deparse(body(f)))
}
saveRDS(definitions,"provenance/runtime-functions.rds")
packages <- c("divextra","diversitree","ape","callr","yaml")
versions <- data.frame(package=packages,version=NA_character_)
for(i in seq_along(packages)) versions$version[i] <- as.character(packageVersion(packages[i]))
write.csv(versions,"provenance/package-versions.csv",row.names=FALSE)
writeLines(capture.output(sessionInfo()),"provenance/local-session.txt")
cat("Validated 15 groups, 60 likelihood comparisons and nine optimization starts.\n")
