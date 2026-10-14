# Run from the DivExtra repository root. Independent fixed-time recovery pilot.
# Rscript inst/benchmarks/archive/pilots/geosse-parameter-recovery/R/06-fixed-time-pilot.R
library(divextra)
library(ggplot2)
source('inst/benchmarks/archive/pilots/geosse-parameter-recovery/R/benchmark-functions.R')
out <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/fixed-time-pilot'
dir.create(out, recursive = TRUE, showWarnings = FALSE)
scenarios <- read_benchmark_scenarios(FALSE)[1:5, ]
scenarios$scenario_number <- 1:5

# Simulate at a fixed time; return a reason for every excluded attempt. ----
simulate.fixed <- function(scenario, old.time, recent.time, seed) {
  set.seed(seed)
  root <- sample(1:3, 1)
  pars <- make_simulation_parameters(scenario)
  tb <- make.tree.classe.td(pars$classe, k=3, max.taxa=3000,
    max.t1=old.time, max.t2=old.time+recent.time, x0=root,
    single.lineage=TRUE)
  if (is.null(tb)) return(list(status='no_tree'))
  if (abs(attr(tb,'t.epoch1')-old.time)>1e-8 ||
      abs(attr(tb,'t.total')-old.time-recent.time)>1e-8)
    return(list(status='taxon_limit'))
  tree <- tryCatch(table2tree(tb), error=function(e) NULL)
  if (is.null(tree) || ape::Ntip(tree)<2) return(list(status='no_tree'))
  depth <- ape::node.depth.edgelength(tree)
  height <- max(depth)
  cut <- height-recent.time
  if (cut<=0) return(list(status='no_old_crown'))
  a <- depth[tree$edge[,1]]
  b <- depth[tree$edge[,2]]
  old.length <- sum(pmax(0,pmin(b,cut)-a))
  diagnostics <- data.frame(N=ape::Ntip(tree), boundary.lineages=sum(a<cut & b>=cut),
    old.nodes=sum(depth[(ape::Ntip(tree)+1):length(depth)]<cut),
    old.length=old.length, recent.length=sum(tree$edge.length)-old.length,
    crown.age=height, boundary.fraction=recent.time/height)
  list(status='accepted',tree=tree,table=tb,diagnostics=diagnostics,
    root=root,pars=pars,old.time=old.time,recent.time=recent.time,seed=seed)
}

# Choose durations using independent calibration trees, before fitting. ----
grid <- expand.grid(old.time=c(15,20,25,30,35),recent.time=c(3,6))
calibration <- list()
design <- scenarios
for (s in 1:5) {
  for (g in seq_len(nrow(grid))) {
    info <- list()
    status <- character()
    for (r in 1:20) {
      z <- simulate.fixed(scenarios[s,],grid$old.time[g],grid$recent.time[g],
        8100000+s*10000+g*100+r)
      status[r] <- z$status
      if (z$status=='accepted') info[[length(info)+1]] <- z$diagnostics
    }
    d <- do.call(rbind,info)
    row <- data.frame(scenario=scenarios$scenario_id[s],old.time=grid$old.time[g],
      recent.time=grid$recent.time[g],accepted=sum(status=='accepted'),
      capped=sum(status=='taxon_limit'),no.tree=sum(status=='no_tree'),
      no.old=sum(status=='no_old_crown'),median.N=NA_real_,median.cross=NA_real_,score=Inf)
    if (!is.null(d) && nrow(d)>0) {
      row$median.N <- median(d$N)
      row$median.cross <- median(d$boundary.lineages)
      if (row$accepted>=8 && row$capped==0)
        row$score <- abs(log(row$median.cross/30))+0.5*abs(log(row$median.N/200))
    }
    calibration[[length(calibration)+1]] <- row
  }
  tab <- do.call(rbind,calibration)
  candidates <- tab[tab$scenario==scenarios$scenario_id[s],]
  best <- candidates[which.min(candidates$score),]
  stopifnot(is.finite(best$score))
  design$old.time[s] <- best$old.time
  design$recent.time[s] <- best$recent.time
  message('Design ',design$scenario_id[s],': old=',best$old.time,
    ', recent=',best$recent.time,', median N=',best$median.N,
    ', boundary lineages=',best$median.cross)
}
write.csv(do.call(rbind,calibration),file.path(out,'calibration.csv'),row.names=FALSE)
write.csv(design,file.path(out,'design.csv'),row.names=FALSE)

# Generate five independent retained trees per scenario; log all attempts. ----
datasets <- list()
attempts <- list()
for (s in 1:5) {
  accepted <- 0
  for (attempt in 1:500) {
    z <- simulate.fixed(design[s,],design$old.time[s],design$recent.time[s],
      9200000+s*10000+attempt)
    attempts[[length(attempts)+1]] <- data.frame(scenario=design$scenario_id[s],
      attempt=attempt,status=z$status)
    if (z$status!='accepted') next
    accepted <- accepted+1
    z$id <- paste0(design$scenario_id[s],'-',accepted)
    z$scenario <- design$scenario_id[s]
    z$replicate <- accepted
    datasets[[length(datasets)+1]] <- z
    if (accepted==5) break
  }
  stopifnot(accepted==5)
}
saveRDS(datasets,file.path(out,'datasets.rds'))
write.csv(do.call(rbind,attempts),file.path(out,'simulation-attempts.csv'),row.names=FALSE)

# One optimizer per process: three distinct, saved starting vectors. ----
fit.start <- function(dataset,start) {
  library(divextra)
  cfg <- read_yaml_pars_td('inst/benchmarks/archive/pilots/geosse-parameter-recovery/config/model-episodic.yml')
  cfg$epoch.times <- dataset$recent.time
  lik <- make.classe.td(dataset$tree,dataset$tree$tip.state,k=3,n.epoch=2,
    sampling.f=c(1,1,1),strict=FALSE,control=list(backend='gslode'))
  lik <- diversitree::constrain(lik,formulae=make_constraints_sse_td(cfg))
  initial <- init.pars.classe_td(lik,dataset$tree,k=3,n.epoch=2,eps=c(.1,.4,.8)[start])
  set.seed(dataset$seed+start*100000)
  initial <- initial*exp(runif(length(initial),-log(2),log(2)))
  warnings <- character()
  fit <- withCallingHandlers(diversitree::find.mle(lik,initial,
    root=diversitree::ROOT.OBS,condition.surv=TRUE,keep.func=FALSE,
    control=list(maxit=10000)),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w)); invokeRestart('muffleWarning')
    })
  list(fit=fit,initial=initial,warnings=unique(warnings))
}
jobs <- expand.grid(dataset=seq_along(datasets),start=1:3)
fits <- list()
for (first in seq(1,nrow(jobs),by=3)) {
  indices <- first:min(first+2,nrow(jobs))
  processes <- list()
  for (j in indices) {
    processes[[as.character(j)]] <- callr::r_bg(fit.start,
      args=list(dataset=datasets[[jobs$dataset[j]]],start=jobs$start[j]))
  }
  started <- Sys.time()
  repeat {
    alive <- FALSE
    for (j in indices) {
      p <- processes[[as.character(j)]]
      if (p$is_alive()) {
        if (as.numeric(difftime(Sys.time(),started,units='secs'))>90) p$kill()
        else alive <- TRUE
      }
    }
    if (!alive) break
    Sys.sleep(.2)
  }
  for (j in indices) {
    result <- tryCatch(processes[[as.character(j)]]$get_result(),error=function(e) list(error=conditionMessage(e)))
    result$id <- datasets[[jobs$dataset[j]]]$id
    result$start <- jobs$start[j]
    fits[[j]] <- result
    message('Finished ',result$id,' start ',result$start,
      if (is.null(result$error)) paste0(' convergence=',result$fit$convergence) else ' failed/timeout')
  }
  saveRDS(fits,file.path(out,'fits.rds'))
}

# Recovery and convergence diagnostics. ----
estimates <- list()
checks <- list()
starts <- list()
for (i in seq_along(datasets)) {
  z <- datasets[[i]]
  ff <- fits[jobs$dataset==i]
  ll <- rep(NA_real_,3)
  converged <- rep(FALSE,3)
  for (j in 1:3) {
    f <- ff[[j]]
    if (is.null(f$error)) {
      ll[j] <- f$fit$lnLik
      converged[j] <- isTRUE(f$fit$convergence==0) && is.finite(ll[j])
    }
    starts[[length(starts)+1]] <- data.frame(id=z$id,start=j,
      converged=converged[j],lnLik=ll[j],error=if(is.null(f$error)) '' else f$error)
  }
  valid <- which(converged)
  spread <- if(length(valid)>=2) diff(range(ll[valid])) else NA_real_
  parameter.spread <- NA_real_
  if (length(valid)>=2) {
    pp <- matrix(NA_real_,length(valid),6)
    for (v in seq_along(valid)) pp[v,] <- ff[[valid[v]]]$fit$par
    parameter.spread <- max(apply(pp,2,function(x) diff(range(x))/max(abs(x),1e-6)))
  }
  checks[[i]] <- cbind(data.frame(id=z$id,scenario=z$scenario,
    converged.starts=length(valid),lnLik.spread=spread,
    parameter.relative.spread=parameter.spread,
    all.starts.agree=length(valid)==3 && spread<1e-4),z$diagnostics)
  if (!length(valid)) next
  best <- valid[which.max(ll[valid])]
  pars <- extract_epoch_estimates(list(best=ff[[best]]$fit))
  for (name in names(pars)) {
    bits <- strsplit(name,'_')[[1]]
    truth <- z$pars$geosse[bits[2],switch(bits[1],s='sA',d='dA',x='xA')]
    estimates[[length(estimates)+1]] <- data.frame(id=z$id,scenario=z$scenario,
      parameter=name,truth=truth,estimate=unname(pars[name]),N=z$diagnostics$N,
      boundary.lineages=z$diagnostics$boundary.lineages)
  }
}
checks <- do.call(rbind,checks)
estimates <- do.call(rbind,estimates)
write.csv(checks,file.path(out,'convergence.csv'),row.names=FALSE)
write.csv(do.call(rbind,starts),file.path(out,'starts.csv'),row.names=FALSE)
write.csv(estimates,file.path(out,'estimates.csv'),row.names=FALSE)
estimates$ratio <- estimates$estimate/estimates$truth
estimates$parameter <- factor(estimates$parameter,
  levels=c('s_recent','d_recent','x_recent','s_old','d_old','x_old'))
p <- ggplot(estimates,aes(scenario,ratio,color=scenario))+
  geom_hline(yintercept=1,linetype=2,color='grey50')+
  geom_point(position=position_jitter(width=.12,height=0,seed=1))+
  facet_wrap(~parameter,ncol=3)+scale_y_log10()+theme_classic()+
  theme(legend.position='none')+labs(x='Scenario',y='Estimate / truth (log scale)',
    caption='Five fixed-time datasets per scenario. Best converged start; failed starts remain in diagnostics.')
ggsave(file.path(out,'recovery.png'),p,width=9,height=6,dpi=180)
ggsave(file.path(out,'recovery.pdf'),p,width=9,height=6)
p2 <- ggplot(checks,aes(N,boundary.lineages,color=scenario))+
  geom_point(size=3)+scale_x_log10()+scale_y_log10()+theme_classic()+
  labs(x='Extant tips',y='Reconstructed lineages at epoch boundary')
ggsave(file.path(out,'epoch-information.png'),p2,width=7,height=5,dpi=180)
ggsave(file.path(out,'epoch-information.pdf'),p2,width=7,height=5)
print(checks)
