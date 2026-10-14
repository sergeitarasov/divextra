# Profile P05 rates, re-optimizing all five nuisance rates at each grid point.
library(divextra)
library(ggplot2)
base <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/fixed-time-pilot'
out <- file.path(base,'profiles')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
datasets <- readRDS(file.path(base,'datasets.rds'))
fits <- readRDS(file.path(base,'fits.rds'))
ids <- c('P05-1','P05-3','P05-4')
tasks <- list()
for (id in ids) {
  z <- datasets[[which(vapply(datasets,function(x) x$id==id,logical(1)))]]
  ff <- fits[vapply(fits,function(x) x$id==id,logical(1))]
  best <- ff[[which.max(vapply(ff,function(x) x$fit$lnLik,numeric(1)))]]$fit
  for (rate in c('q13.2','mu1.2')) {
    grid <- if(rate=='q13.2') c(1e-8,.003,.01,.03,.1,.3,1,3,10) else c(1e-8,.002,.006,.02,.06,.2,.6,2)
    grid <- sort(unique(c(grid,unname(best$par[rate]))))
    for (value in grid) tasks[[length(tasks)+1]] <- list(z=z,rate=rate,value=value,best=best)
  }
}
profile.one <- function(task) {
  library(divextra)
  cfg <- read_yaml_pars_td('inst/benchmarks/archive/pilots/geosse-parameter-recovery/config/model-episodic.yml')
  cfg$epoch.times <- task$z$recent.time
  build <- function(backend,tol) {
    lik <- make.classe.td(task$z$tree,task$z$tree$tip.state,k=3,n.epoch=2,
      sampling.f=c(1,1,1),strict=FALSE,control=list(backend=backend,tol=tol))
    diversitree::constrain(lik,formulae=make_constraints_sse_td(cfg))
  }
  lik <- build('gslode',1e-10)
  full <- task$best$par
  free <- names(full)!=task$rate
  truth <- c(lambda111.1=.30,mu1.1=.15,q13.1=.20,
    lambda111.2=.10,mu1.2=.02,q13.2=.03)
  objective <- function(logp) {
    p <- full
    p[free] <- exp(logp)
    p[task$rate] <- task$value
    v <- tryCatch(suppressWarnings(lik(p,root=diversitree::ROOT.OBS,condition.surv=TRUE)),error=function(e) -Inf)
    if (!is.finite(v)) return(1e50)
    -v
  }
  answers <- list()
  for (s in 1:2) {
    initial <- if(s==1) full[free] else truth[names(full)[free]]
    ans <- nlminb(log(pmax(initial,1e-9)),objective,lower=log(1e-10),upper=log(100),
      control=list(iter.max=500,eval.max=1500,rel.tol=1e-9))
    answers[[s]] <- ans
  }
  chosen <- which.min(vapply(answers,function(x) x$objective,numeric(1)))
  opt <- answers[[chosen]]
  full[free] <- exp(opt$par)
  full[task$rate] <- task$value
  vals <- numeric()
  for (backend in c('gslode','deSolve')) {
    for (tol in c(1e-8,1e-10,1e-12)) {
      f <- build(backend,tol)
      vals[paste(backend,tol)] <- tryCatch(suppressWarnings(f(full,
        root=diversitree::ROOT.OBS,condition.surv=TRUE)),error=function(e) NA_real_)
    }
  }
  list(id=task$z$id,rate=task$rate,value=task$value,lnLik=-opt$objective,
    convergence=opt$convergence,parameters=full,optimizers=answers,numerical=vals,
    bound.hit=any(full[free]>99))
}
results <- list()
for (first in seq(1,length(tasks),by=3)) {
  indices <- first:min(first+2,length(tasks))
  processes <- list()
  for (j in indices) processes[[as.character(j)]] <- callr::r_bg(profile.one,list(task=tasks[[j]]))
  time <- Sys.time()
  repeat {
    alive <- FALSE
    for (j in indices) {
      p <- processes[[as.character(j)]]
      if(p$is_alive()) {
        if(as.numeric(difftime(Sys.time(),time,units='secs'))>120) p$kill() else alive <- TRUE
      }
    }
    if(!alive) break
    Sys.sleep(.2)
  }
  for(j in indices) {
    results[[j]] <- tryCatch(processes[[as.character(j)]]$get_result(),
      error=function(e) list(id=tasks[[j]]$z$id,rate=tasks[[j]]$rate,value=tasks[[j]]$value,error=conditionMessage(e)))
    message(j,'/',length(tasks),' ',tasks[[j]]$z$id,' ',tasks[[j]]$rate,'=',tasks[[j]]$value,
      if(is.null(results[[j]]$error)) paste(' logLik',round(results[[j]]$lnLik,3)) else ' failed')
  }
  saveRDS(results,file.path(out,'profiles.rds'))
}
rows <- list()
for(r in results) {
  if(!is.null(r$error)) next
  rows[[length(rows)+1]] <- data.frame(id=r$id,rate=r$rate,value=r$value,lnLik=r$lnLik,
    convergence=r$convergence,bound.hit=r$bound.hit,
    numerical.range=diff(range(r$numerical,na.rm=TRUE)),
    numerical.finite=sum(is.finite(r$numerical)))
}
d <- do.call(rbind,rows)
d$loss <- NA_real_
for(id in ids) {
  keep <- d$id==id
  old <- fits[vapply(fits,function(x) x$id==id,logical(1))]
  reference <- max(c(d$lnLik[keep],vapply(old,function(x) x$fit$lnLik,numeric(1))))
  d$loss[keep] <- reference-d$lnLik[keep]
}
write.csv(d,file.path(out,'profiles.csv'),row.names=FALSE)
truth <- expand.grid(id=ids,rate=c('q13.2','mu1.2'))
truth$value <- ifelse(truth$rate=='q13.2',.03,.02)
p <- ggplot(d,aes(value,loss)) + geom_line()+geom_point(aes(shape=factor(convergence)))+
  geom_vline(data=truth,aes(xintercept=value),color='#0072B2',linetype=2)+
  geom_hline(yintercept=1.92,linetype=3,color='grey40')+
  facet_grid(id~rate,scales='free_x')+scale_x_log10()+theme_classic()+
  labs(x='Fixed old-epoch rate (log scale)',y='Log-likelihood loss',shape='Optimizer code',
    caption='q13.2 = dispersal; mu1.2 = contraction/extinction. Blue: truth. 1.92 is a reference, not calibrated coverage.')
ggsave(file.path(out,'profiles.png'),p,width=9,height=8,dpi=180)
ggsave(file.path(out,'profiles.pdf'),p,width=9,height=8)
print(d)
