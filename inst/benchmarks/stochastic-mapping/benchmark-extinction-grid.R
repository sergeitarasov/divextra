# Numerical extinction-cache benchmark ----
# Pure-death stretched-grid prototype only: NOT a general ClaSSE method.
# Production code is unchanged. Load the development package before sourcing
# this file when testing a checkout; otherwise use the installed package.
if (!'divextra' %in% loadedNamespaces()) library(divextra)
cache.fun <- divextra:::.classe_td_extinction_cache
schedule.fun <- divextra:::.classe_td_schedule
results <- list()
for (mu in c(.7,100,10000)) {
  sch <- schedule.fun(c(lambda111.1=0,mu1.1=mu),1,'A')
  # Fixed grids are benchmarked at bounded age to avoid multi-million-row caches.
  T <- min(3,2.5/mu)
  for (h in unique(c(.01,.1/mu,.01/mu))) {
    elapsed <- system.time(cc <- tryCatch(cache.fun(sch,T,list(max.step=h)),error=function(e)e))[['elapsed']]
    if (inherits(cc,'error')) {
      results[[length(results)+1]] <- data.frame(mu=mu,T=T,h=h,seconds=elapsed,points=NA,E.error=NA,CDF.error=NA,status=conditionMessage(cc))
    } else {
      tt <- sort(unique(c(seq(0,T,length.out=20001),exp(seq(log(T*1e-9),log(T),length.out=10000)))))
      tt <- pmin(T,tt)
      ex <- -expm1(-mu*tt)
      ap <- vapply(tt,cc$Ei,numeric(1),state=1)
      cdf.ex <- exp(-mu*(T-tt))*ex/(-expm1(-mu*T))
      cdf.ap <- exp(-mu*(T-tt))*ap/cc$Ei(T,1)
      results[[length(results)+1]] <- data.frame(mu=mu,T=T,h=h,seconds=elapsed,points=length(cc$times),E.error=max(abs(ap-ex)),CDF.error=max(abs(cdf.ex-cdf.ap)),status='ok')
    }
  }
}
print(do.call(rbind,results),digits=7,row.names=FALSE)

# Pure-death proof-of-concept: stretched grid, not a general ClaSSE algorithm.
# Substitute only local clone's mesh; production function untouched.
stretched <- cache.fun
bb <- paste(deparse(body(stretched),width.cutoff=500),collapse='\n')
needle <- 'seq(lo, hi, length.out = ceiling(span/control$max.step) + 1L)'
stopifnot(grepl(needle,bb,fixed=TRUE))
bb <- gsub(needle,'lo + stretch.mesh(span, max(totals[, epoch]))',bb,fixed=TRUE)
body(stretched) <- parse(text=bb)[[1]]
ee <- new.env(parent=environment(cache.fun))
ee$stretch.mesh <- function(span,r) {
  eta <- .003
  pmin(span,expm1(seq(0,log1p(r*span),length.out=ceiling(log1p(r*span)/log1p(eta))+1))/r)
}
environment(stretched)<-ee
out <- list()
for(mu in c(.7,100,10000)) {
  sch <- schedule.fun(c(lambda111.1=0,mu1.1=mu),1,'A')
  elapsed<-system.time(cc<-tryCatch(stretched(sch,3,list(rtol=1e-12)),error=function(e)e))[['elapsed']]
  if(inherits(cc,'error')) { print(conditionMessage(cc)); next }
  tt<-sort(unique(c(seq(0,3,length.out=10000),exp(seq(log(1e-12),log(3),length.out=20000)))))
  tt<-pmin(3,tt)
  ex<--expm1(-mu*tt); ap<-vapply(tt,cc$Ei,numeric(1),state=1)
  a<-min(3,2.5/mu); uu<-seq(0,a,length.out=10000)
  exact<-exp(-mu*(a-uu))*(-expm1(-mu*uu))/(-expm1(-mu*a))
  approx<-exp(-mu*(a-uu))*vapply(uu,cc$Ei,numeric(1),state=1)/cc$Ei(a,1)
  out[[length(out)+1]]<-data.frame(mu=mu,T=3,seconds=elapsed,points=length(cc$times),E.error=max(abs(ap-ex)),CDF.error=max(abs(exact-approx)))
}
print(do.call(rbind,out),digits=7,row.names=FALSE)
