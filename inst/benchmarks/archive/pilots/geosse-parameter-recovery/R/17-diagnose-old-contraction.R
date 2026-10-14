# Diagnose old contraction without changing inference or previous results. ----
library(divextra)
library(ggplot2)
base <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/epoch-grid-pilot'
out <- file.path(base,'old-contraction-diagnostic')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
datasets <- readRDS(file.path(base,'datasets.rds'))
best <- readRDS(file.path(base,'best-fits.rds'))
design <- read.csv(file.path(base,'design.csv'))
estimates <- read.csv(file.path(base,'estimates.csv'))

# Full-history oracle: x acts once in A/B and twice in AB. ----
oracle <- list()
selected <- list()
for(s in design$scenario_id[design$rate.scenario %in% c('all_rates_change','contraction_only')]) {
  e <- estimates[estimates$scenario==s & estimates$label=='x_old',]
  target <- e$id[order(e$N)][3]
  for(z in datasets) if(z$scenario==s) {
    row <- design[design$scenario_id==s,]
    tb <- z$table
    hist <- attr(tb,'hist')
    exposure <- c(single=0,AB=0)
    for(i in seq_len(nrow(tb))) {
      h <- hist[hist$idx==tb$idx[i],,drop=FALSE]
      h <- h[order(h$t),,drop=FALSE]
      times <- c(tb$start[i],h$t,tb$start[i]+tb$len[i])
      states <- if(nrow(h)) c(h$from[1],h$to) else tb$state[i]
      lengths <- pmax(0,pmin(times[-1],row$old.time)-times[-length(times)])
      exposure['single'] <- exposure['single']+sum(lengths[states!=3])
      exposure['AB'] <- exposure['AB']+sum(lengths[states==3])
    }
    deaths <- z$diagnostics$old.sim.extinction
    contractions <- z$diagnostics$old.sim.contraction
    oracle[[length(oracle)+1]] <- data.frame(id=z$id,scenario=s,N=z$diagnostics$N,
      deaths=deaths,contractions=contractions,single.exposure=exposure['single'],
      AB.exposure=exposure['AB'],oracle.x=(deaths+contractions)/(exposure['single']+2*exposure['AB']),
      reconstructed.x=best[[paste(z$id,'linked_classe')]]$par['x_old'])
    if(z$id==target) selected[[length(selected)+1]] <- z
  }
}
oracle <- do.call(rbind,oracle)
write.csv(oracle,file.path(out,'full-history-diagnostic.csv'),row.names=FALSE)

# Profile x_old, optimizing all other independent parameters. ----
profile.one <- function(z,best,value) {
  library(divextra)
  cfg <- read_yaml_pars_td('inst/benchmarks/archive/pilots/geosse-parameter-recovery/config/model-episodic.yml')
  cfg$epoch.times <- z$recent.time
  build <- function(tol) {
    raw <- make.classe.td(z$tree,z$tree$tip.state,k=3,n.epoch=2,strict=FALSE,
      sampling.f=c(1,1,1),control=list(tol=tol))
    diversitree::constrain(raw,formulae=make_constraints_sse_td(cfg))
  }
  raw <- build(1e-10)
  expand <- function(p) {
    full <- z$truth
    for(rate in c('s','d','x')) {
      full[paste0(rate,'_recent')] <- p[paste0(rate,'_recent')]
      source <- if(z$links[rate]) paste0(rate,'_recent') else paste0(rate,'_old')
      full[paste0(rate,'_old')] <- p[source]
    }
    full[c(1,3,2,4,6,5)]
  }
  evaluate <- function(p) raw(expand(p),root=diversitree::ROOT.OBS,condition.surv=TRUE)
  free <- setdiff(names(best$par),'x_old')
  objective <- function(logp) {
    p <- best$par
    p[free] <- exp(logp)
    p['x_old'] <- value
    ll <- tryCatch(suppressWarnings(evaluate(p)),error=function(e) -Inf)
    if(!is.finite(ll)) return(1e50)
    -ll
  }
  answers <- list()
  for(start in 1:2) {
    initial <- if(start==1) best$par[free] else z$truth[free]
    answers[[start]] <- nlminb(log(pmax(initial,1e-10)),objective,
      lower=log(1e-10),upper=log(100),control=list(iter.max=800,eval.max=2500,rel.tol=1e-10))
  }
  winner <- 1
  for(j in 2:length(answers)) if(answers[[j]]$objective<answers[[winner]]$objective) winner <- j
  fit <- answers[[winner]]
  pars <- best$par
  pars[free] <- exp(fit$par)
  pars['x_old'] <- value
  numerical <- numeric()
  for(tol in c(1e-8,1e-10,1e-12)) {
    test <- build(tol)
    numerical <- c(numerical,test(expand(pars),root=diversitree::ROOT.OBS,condition.surv=TRUE))
  }
  truth <- z$truth[names(best$par)]
  truth['x_old'] <- value
  list(id=z$id,value=value,lnLik=-fit$objective,convergence=fit$convergence,
    par=pars,numerical.range=diff(range(numerical)),
    oracle.nuisance.lnLik=evaluate(truth),optimizers=answers)
}
tasks <- list()
for(z in selected) for(value in c(0,.002,.01,.02,.05,.1,.2))
  tasks[[length(tasks)+1]] <- list(z=z,best=best[[paste(z$id,'linked_classe')]],value=value)
results <- vector('list',length(tasks))
if(file.exists(file.path(out,'profiles.rds'))) results <- readRDS(file.path(out,'profiles.rds'))
pending <- which(vapply(results,is.null,logical(1)))
if(length(pending)) for(first in seq(1,length(pending),by=3)) {
  indices <- pending[first:min(first+2,length(pending))]
  workers <- list()
  for(j in indices) workers[[as.character(j)]] <- callr::r_bg(profile.one,tasks[[j]])
  begin <- Sys.time()
  repeat {
    alive <- FALSE
    for(j in indices) if(workers[[as.character(j)]]$is_alive()) {
      if(as.numeric(difftime(Sys.time(),begin,units='secs'))>120) workers[[as.character(j)]]$kill() else alive <- TRUE
    }
    if(!alive) break
    Sys.sleep(.2)
  }
  for(j in indices) {
    results[[j]] <- tryCatch(workers[[as.character(j)]]$get_result(),error=function(e)
      list(id=tasks[[j]]$z$id,value=tasks[[j]]$value,error=conditionMessage(e)))
    message(j,'/',length(tasks),' ',tasks[[j]]$z$id,' x=',tasks[[j]]$value)
  }
  saveRDS(results,file.path(out,'profiles.rds'))
}
rows <- list()
for(r in results) {
  if(!is.null(r$error)) next
  rows[[length(rows)+1]] <- data.frame(id=r$id,x=r$value,lnLik=r$lnLik,
    convergence=r$convergence,numerical.range=r$numerical.range,
    oracle.nuisance.lnLik=r$oracle.nuisance.lnLik)
}
pdata <- do.call(rbind,rows)
for(id in unique(pdata$id)) {
  keep <- pdata$id==id
  ref <- max(c(pdata$lnLik[keep],best[[paste(id,'linked_classe')]]$lnLik))
  pdata$loss[keep] <- ref-pdata$lnLik[keep]
  pdata$oracle.nuisance.loss[keep] <- max(pdata$oracle.nuisance.lnLik[keep])-pdata$oracle.nuisance.lnLik[keep]
}
write.csv(pdata,file.path(out,'profiles.csv'),row.names=FALSE)
p <- ggplot(pdata,aes(x,loss))+geom_line()+geom_point(aes(shape=factor(convergence)))+
  geom_hline(yintercept=1.92,linetype=2)+geom_vline(xintercept=.02,linetype=3)+
  facet_wrap(~id,scales='free_y')+theme_classic()+
  labs(x='Old contraction rate',y='Profile log-likelihood loss',shape='Optimizer code',
    caption='Dotted: truth. Dashed: 1.92 reference only, not calibrated coverage near a boundary. Two nuisance starts per point.')
ggsave(file.path(out,'profiles.png'),p,width=11,height=6,dpi=160)
ggsave(file.path(out,'profiles.pdf'),p,width=11,height=6)
print(pdata[pdata$x %in% c(0,.02),],row.names=FALSE)
print(aggregate(cbind(deaths,contractions,oracle.x,reconstructed.x)~scenario,oracle,median),row.names=FALSE)
