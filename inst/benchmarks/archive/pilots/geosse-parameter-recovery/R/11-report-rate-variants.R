# Report recovery, convergence and the classical GeoSSE equivalence check.
library(divextra)
library(ggplot2)
out <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/rate-variants-pilot'
datasets <- readRDS(file.path(out,'datasets.rds'))
fits <- readRDS(file.path(out,'fits.rds'))
design <- read.csv(file.path(out,'design.csv'))
rows <- list()
for(f in fits) {
  if(is.null(f)) next
  ok <- is.null(f$error) && isTRUE(f$convergence==0) && is.finite(f$lnLik)
  rows[[length(rows)+1]] <- data.frame(id=f$id,scenario=f$scenario,model=f$model,start=f$start,
    converged=ok,lnLik=if(is.null(f$error)) f$lnLik else NA_real_,
    convergence.code=if(is.null(f$error)) f$convergence else NA_integer_,
    message=if(is.null(f$error)) f$message else f$error,
    bound.hit=if(is.null(f$error)) f$bound.hit else NA,
    elapsed=if(is.null(f$error)) f$elapsed else NA_real_,
    error=if(is.null(f$error)) '' else f$error)
}
starts <- do.call(rbind,rows)
write.csv(starts,file.path(out,'starts.csv'),row.names=FALSE)
checks <- list()
estimates <- list()
best.fits <- list()
for(z in datasets) {
  for(model in unique(starts$model[starts$id==z$id])) {
    ss <- starts[starts$id==z$id & starts$model==model,]
    valid <- ss[ss$converged,]
    key <- paste(z$id,model)
    span <- if(nrow(valid)>1) diff(range(valid$lnLik)) else NA_real_
    parameter.spread <- NA_real_
    vectors <- list()
    for(f in fits) if(!is.null(f) && f$id==z$id && f$model==model && f$start %in% valid$start)
      vectors[[length(vectors)+1]] <- f$par
    if(length(vectors)>1) {
      matrix.par <- do.call(rbind,vectors)
      spreads <- numeric(ncol(matrix.par))
      for(column in seq_len(ncol(matrix.par))) spreads[column] <-
        diff(range(matrix.par[,column]))/max(abs(matrix.par[,column]),1e-6)
      parameter.spread <- max(spreads)
    }
    checks[[length(checks)+1]] <- cbind(data.frame(id=z$id,scenario=z$scenario,model=model,
      converged.starts=nrow(valid),lnLik.spread=span,
      parameter.relative.spread=parameter.spread,
      all.agree=nrow(valid)==3 && span<1e-4),z$diagnostics)
    if(!nrow(valid)) next
    winner <- valid$start[which.max(valid$lnLik)]
    for(f in fits) if(!is.null(f) && f$id==z$id && f$model==model && f$start==winner) best.fits[[key]] <- f
    best <- best.fits[[key]]
    if(model=='episodic') for(p in names(z$truth)) estimates[[length(estimates)+1]] <- data.frame(
      id=z$id,scenario=z$scenario,parameter=p,truth=unname(z$truth[p]),
      estimate=unname(best$par[p]),N=z$diagnostics$N,boundary.lineages=z$diagnostics$boundary.lineages)
  }
}
checks <- do.call(rbind,checks)
estimates <- do.call(rbind,estimates)
write.csv(checks,file.path(out,'convergence.csv'),row.names=FALSE)
write.csv(estimates,file.path(out,'estimates.csv'),row.names=FALSE)
saveRDS(best.fits,file.path(out,'best-fits.rds'))

# Compare constant ClaSSE and GeoSSE at exactly matched parameter vectors. ----
audit <- list()
comparison <- list()
for(z in datasets) {
  epi <- best.fits[[paste(z$id,'episodic')]]
  cl <- best.fits[[paste(z$id,'constant_classe')]]
  row <- data.frame(id=z$id,scenario=z$scenario,
    delta.AIC.constant.minus.episodic=NA_real_,nested.likelihood.violation=NA,
    episodic.lnLik.minus.truth=if(is.null(epi)) NA_real_ else epi$lnLik-epi$truth.lnLik)
  if(!is.null(epi) && !is.null(cl)) {
    row$delta.AIC.constant.minus.episodic <- (6-2*cl$lnLik)-(12-2*epi$lnLik)
    row$nested.likelihood.violation <- epi$lnLik<cl$lnLik-1e-4
  }
  comparison[[length(comparison)+1]] <- row
  if(z$scenario!='constant') next
  geo <- best.fits[[paste(z$id,'geosse')]]
  states <- setNames(c(1L,2L,0L)[z$tree$tip.state],names(z$tree$tip.state))
  g <- diversitree::make.geosse(z$tree,states,strict=FALSE,control=list(tol=1e-10))
  cl.likelihood <- make.classe.td(z$tree,z$tree$tip.state,k=3,n.epoch=2,strict=FALSE,control=list(tol=1e-10))
  probes <- list(truth=c(.15,.10,.05),low=c(.05,.03,.01),high=c(.4,.3,.2))
  if(!is.null(cl)) probes$classe.optimum <- unname(cl$par)
  if(!is.null(geo)) probes$geosse.optimum <- unname(geo$par)
  for(probe in names(probes)) {
    p <- probes[[probe]]
    gp <- c(sA=p[1],sB=p[1],sAB=p[1],xA=p[3],xB=p[3],dA=p[2],dB=p[2])
    cp <- pars.geosse2classe(gp)
    tp <- c(t.1=z$recent.time,setNames(cp,paste0(names(cp),'.1')),setNames(cp,paste0(names(cp),'.2')))
    for(survival in c(FALSE,TRUE)) {
      ll.g <- g(gp,root=diversitree::ROOT.OBS,condition.surv=survival)
      ll.c <- cl.likelihood(tp,root=diversitree::ROOT.OBS,condition.surv=survival)
      audit[[length(audit)+1]] <- data.frame(id=z$id,probe=probe,survival=survival,
        geosse=ll.g,classe=ll.c,difference=ll.c-ll.g)
    }
  }
}
comparison <- do.call(rbind,comparison)
audit <- do.call(rbind,audit)
write.csv(comparison,file.path(out,'model-comparison.csv'),row.names=FALSE)
write.csv(audit,file.path(out,'geosse-likelihood-audit.csv'),row.names=FALSE)

control <- list()
for(z in datasets) if(z$scenario=='constant') {
  cl <- best.fits[[paste(z$id,'constant_classe')]]
  ge <- best.fits[[paste(z$id,'geosse')]]
  if(is.null(cl) || is.null(ge)) next
  for(p in c('s','d','x')) control[[length(control)+1]] <- data.frame(id=z$id,parameter=p,
    truth=switch(p,s=.15,d=.10,x=.05),classe=unname(cl$par[p]),geosse=unname(ge$par[p]),
    lnLik.difference=cl$lnLik-ge$lnLik)
}
control <- do.call(rbind,control)
write.csv(control,file.path(out,'geosse-estimate-comparison.csv'),row.names=FALSE)

recovery <- list()
for(s in design$scenario_id) for(p in unique(estimates$parameter)) {
  d <- estimates[estimates$scenario==s & estimates$parameter==p,]
  if(!nrow(d)) next
  recovery[[length(recovery)+1]] <- data.frame(scenario=s,parameter=p,n=nrow(d),truth=d$truth[1],
    median=median(d$estimate),mean=mean(d$estimate),
    within.twofold=sum(d$estimate/d$truth>=.5 & d$estimate/d$truth<=2),
    median.absolute.log10.error=median(abs(log10(d$estimate/d$truth))))
}
recovery <- do.call(rbind,recovery)
write.csv(recovery,file.path(out,'recovery-summary.csv'),row.names=FALSE)

estimates$ratio <- estimates$estimate/estimates$truth
estimates$scenario <- factor(estimates$scenario,levels=design$scenario_id)
estimates$parameter <- factor(estimates$parameter,levels=c('s_recent','d_recent','x_recent','s_old','d_old','x_old'))
p <- ggplot(estimates,aes(scenario,ratio,color=scenario))+geom_hline(yintercept=1,linetype=2)+
  geom_point(position=position_jitter(width=.12,height=0,seed=1),size=2)+
  facet_wrap(~parameter,ncol=3)+scale_y_log10()+theme_classic()+
  theme(legend.position='none',axis.text.x=element_text(angle=35,hjust=1))+
  labs(x=NULL,y='Estimate / truth (log scale)',caption='Five trees per scenario. Best converged episodic fit. Individual estimates; no coverage claim.')
ggsave(file.path(out,'recovery.png'),p,width=11,height=7,dpi=180)
ggsave(file.path(out,'recovery.pdf'),p,width=11,height=7)
p2 <- ggplot(control,aes(geosse,classe))+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_point(size=2)+facet_wrap(~parameter,scales='free')+theme_classic()+
  labs(x='Classical GeoSSE estimate',y='Constant ClaSSE estimate',caption='Identical trees, state mapping, root treatment and survival conditioning.')
ggsave(file.path(out,'geosse-comparison.png'),p2,width=9,height=3.5,dpi=180)
ggsave(file.path(out,'geosse-comparison.pdf'),p2,width=9,height=3.5)
starts$loss <- NA_real_
for(id in unique(starts$id)) for(model in unique(starts$model[starts$id==id])) {
  keep <- starts$id==id & starts$model==model
  if(any(is.finite(starts$lnLik[keep]))) starts$loss[keep] <- max(starts$lnLik[keep],na.rm=TRUE)-starts$lnLik[keep]
}
p3 <- ggplot(starts,aes(id,factor(start),fill=log10(1+loss)))+geom_tile()+
  facet_wrap(~model,ncol=1)+scale_fill_viridis_c(na.value='grey80')+theme_classic()+
  theme(axis.text.x=element_text(angle=90,vjust=.5,hjust=1))+
  labs(x=NULL,y='Start',fill='log10(1 + loss)',caption='Grey or absent entries have no likelihood. GeoSSE was fitted only to constant-rate control trees.')
ggsave(file.path(out,'convergence.png'),p3,width=12,height=6,dpi=180)
ggsave(file.path(out,'convergence.pdf'),p3,width=12,height=6)
summary <- list()
for(s in design$scenario_id) {
  cc <- checks[checks$scenario==s & checks$model=='episodic',]
  aa <- comparison[comparison$scenario==s,]
  summary[[s]] <- data.frame(scenario=s,median.N=median(cc$N),median.boundary.lineages=median(cc$boundary.lineages),
    episodic.converged.starts=sum(cc$converged.starts),all.starts.agree=sum(cc$all.agree),
    episodic.AIC.support=sum(aa$delta.AIC.constant.minus.episodic>2,na.rm=TRUE))
}
summary <- do.call(rbind,summary)
write.csv(summary,file.path(out,'scenario-summary.csv'),row.names=FALSE)
lines <- c('# Rate variants and classical GeoSSE control','',
  'Five new fixed-time trees per scenario. Each episodic and constant ClaSSE model has three starts.',
  'The five constant control trees are additionally fitted using classical diversitree::make.geosse with sA=sB=sAB, dA=dB, xA=xB.',
  'Parameters use ROOT.OBS, condition.surv=TRUE, complete sampling and tolerance 1e-10.',
  'ClaSSE state order is A, B, AB; GeoSSE uses AB=0, A=1, B=2.',
  'All models use the same nlminb optimizer, log rates, independently perturbed starting vectors and rate bounds [1e-10,100].',
  'Three processes run concurrently; each start has a 120-second wall-time limit. Failed starts remain in starts.csv.',
  'The P05-family durations are held fixed at old=30 and recent=6; constant-control durations are 25 and 6.',
  'No selection on tip count or realized old-epoch information. Trees lacking a reconstructed crown across both epochs are excluded; all exclusions are logged.',
  '', '## Design', '', '```',capture.output(print(design,row.names=FALSE)),'```',
  '', '## Convergence and model comparison', '',
  paste('Converged starts:',sum(starts$converged),'of',nrow(starts)),
  paste('Returned fits at the upper rate bound:',sum(starts$bound.hit,na.rm=TRUE)),
  '', 'Incomplete or nonconverged starts:', '```',
  capture.output(print(starts[!starts$converged,c('id','model','start','convergence.code','message')],row.names=FALSE)),'```',
  paste('Nested likelihood violations (episodic below constant by >0.0001):',sum(comparison$nested.likelihood.violation,na.rm=TRUE)),
  paste('Best episodic fits below the generating-parameter likelihood by >0.0001:',
    sum(comparison$episodic.lnLik.minus.truth< -1e-4,na.rm=TRUE)),
  '', '```',capture.output(print(summary,row.names=FALSE)),'```',
  'Agreement requires all three starts to converge and a log-likelihood range <0.0001. This does not imply accurate parameter recovery.',
  'convergence.csv also records maximum relative parameter spread across converged starts, with a 1e-6 denominator floor for rates near zero.',
  'Episodic AIC support counts datasets where constant AIC minus episodic AIC exceeds 2.',
  '', '## Classical GeoSSE equivalence', '',
  paste('Maximum absolute matched-parameter log-likelihood difference:',format(max(abs(audit$difference)),digits=6)),
  paste('Maximum best-fit log-likelihood difference:',format(max(abs(control$lnLik.difference)),digits=6)),
  paste('Maximum absolute parameter-estimate difference:',format(max(abs(control$classe-control$geosse)),digits=6)),
  'Matched parameter tests include true, low, high and both fitted parameter vectors, with and without survival conditioning.',
  'The GeoSSE and constrained ClaSSE fits describe the same model; they are not counted as separate alternatives in AIC selection.',
  '', '## Recovery', '', '```',capture.output(print(recovery,row.names=FALSE)),'```',
  '', '## Interpretation of this pilot', '',
  'The constant-rate implementation agrees closely with classical GeoSSE, both at matched parameter vectors and after optimization. This is strong evidence for the constant-rate special case, not proof of parameter identifiability or of every epoch-boundary calculation.',
  'Recent speciation and dispersal are recovered much more consistently than old dispersal and contraction. Several old-rate estimates approach zero; one reaches the upper bound. Optimizer convergence alone therefore cannot establish reliable recovery.',
  'Reducing recent speciation did not consistently resolve old-rate recovery. Increasing old dispersal produced the strongest recent-rate recovery, but also much larger trees (median 454 tips versus 75 in the baseline), so this cannot be attributed to the rate change alone.',
  'Increasing old contraction did not resolve its recovery. Even in the higher-old-dispersal scenario, none of the five old-contraction estimates was within a factor of two of truth.',
  'For the five constant controls, AIC preferred the constant model over the six-parameter episodic model. Five controls are too few to estimate a false-positive rate precisely.',
  'For the next benchmark, retain both difficult and informative scenarios, increase replication, report recovery against realized tree size and old-epoch branch exposure, and extend/restart failed or disagreeing optimization runs. Do not treat these pilot results as confidence-interval coverage validation.',
  '', '### Constant-model estimates (not the overparameterized episodic fits)', '',
  '```',capture.output(print(control,row.names=FALSE)),'```',
  '', '## Figures', '',
  '- [Parameter recovery](recovery.pdf): estimate divided by truth; the dashed line is exact recovery. Each point is one tree.',
  '- [Classical GeoSSE comparison](geosse-comparison.pdf): matched constant-model parameter estimates; the dashed line is equality, not the generating truth.',
  '- [Starting-point agreement](convergence.pdf): likelihood loss relative to the best returned start for each tree/model. Consult starts.csv for optimizer status; a finite likelihood is not necessarily converged.',
  '', '## Limits', '',
  'Five replicates per scenario are diagnostic only. Unequal realized tree sizes and old-epoch histories limit comparisons.',
  'These are independent trees across scenarios, not paired common-history simulations. No confidence-interval coverage claim is made.',
  'Changing old dispersal or contraction changes diversification and state occupancy as well as the number of transitions.',
  'Slower recent speciation does not guarantee better recovery. Do not choose a final benchmark only because a pilot gives favorable estimates.',
  'Only the constant control is checked against classical GeoSSE; this does not independently validate every episodic likelihood component.',
  '', '## Reproduction', '',
  'Run R/10-rate-variants-pilot.R and then R/11-report-rate-variants.R from the DivExtra repository root.',
  'The runner reuses saved fits. To change the design use a separate results directory; do not reuse an old fits.rds.',
  '', '```',capture.output(sessionInfo()),'```')
writeLines(lines,file.path(out,'report.md'))
print(summary,row.names=FALSE)
print(control,row.names=FALSE)
