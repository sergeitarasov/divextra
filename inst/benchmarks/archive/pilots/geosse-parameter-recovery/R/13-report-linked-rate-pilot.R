# Report recovery, convergence and the classical GeoSSE equivalence check.
library(divextra)
library(ggplot2)
out <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/linked-rate-pilot'
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
    if(model=='linked_classe') for(p in names(best$par)) estimates[[length(estimates)+1]] <- data.frame(
      id=z$id,scenario=z$scenario,parameter=p,truth=unname(z$truth[p]),
      estimate=unname(best$par[p]),N=z$diagnostics$N,boundary.lineages=z$diagnostics$boundary.lineages,
      old.length=z$diagnostics$old.length,truth.lnLik=best$truth.lnLik,lnLik=best$lnLik)
  }
}
checks <- do.call(rbind,checks)
estimates <- do.call(rbind,estimates)
write.csv(checks,file.path(out,'convergence.csv'),row.names=FALSE)
write.csv(estimates,file.path(out,'estimates.csv'),row.names=FALSE)
saveRDS(best.fits,file.path(out,'best-fits.rds'))


# Recovery of independent parameters (shared rates are counted only once). ----
estimates$rate <- sub('_.*','',estimates$parameter)
estimates$epoch <- sub('.*_','',estimates$parameter)
for(i in seq_len(nrow(estimates))) {
  row <- design[design$scenario_id==estimates$scenario[i],]
  if(row[[paste0('link_',estimates$rate[i])]]) estimates$epoch[i] <- 'shared'
}
estimates$label <- paste(estimates$rate,estimates$epoch,sep='_')
estimates$ratio <- estimates$estimate/estimates$truth
estimates$absolute.log10.error <- abs(log10(estimates$ratio))
write.csv(estimates,file.path(out,'estimates.csv'),row.names=FALSE)
recovery <- list()
for(s in design$scenario_id) for(p in unique(estimates$label)) {
  d <- estimates[estimates$scenario==s & estimates$label==p,]
  if(!nrow(d)) next
  recovery[[length(recovery)+1]] <- data.frame(scenario=s,parameter=p,n=nrow(d),
    truth=d$truth[1],median=median(d$estimate),mean=mean(d$estimate),
    within.twofold=sum(d$ratio>=.5 & d$ratio<=2))
}
recovery <- do.call(rbind,recovery)
write.csv(recovery,file.path(out,'recovery-summary.csv'),row.names=FALSE)

# Classical GeoSSE: same control datasets and same parameterization. ----
control <- list()
audit <- list()
for(z in datasets) if(z$scenario=='constant') {
  cl <- best.fits[[paste(z$id,'linked_classe')]]
  ge <- best.fits[[paste(z$id,'geosse')]]
  if(!is.null(cl) && !is.null(ge)) for(p in c('s','d','x')) {
    control[[length(control)+1]] <- data.frame(id=z$id,parameter=p,
      truth=z$truth[paste0(p,'_recent')],classe=cl$par[paste0(p,'_recent')],
      geosse=ge$par[p],lnLik.difference=cl$lnLik-ge$lnLik)
  }
  states <- setNames(c(1L,2L,0L)[z$tree$tip.state],names(z$tree$tip.state))
  g <- diversitree::make.geosse(z$tree,states,strict=FALSE,control=list(tol=1e-10))
  cl.likelihood <- make.classe.td(z$tree,z$tree$tip.state,k=3,n.epoch=2,strict=FALSE,control=list(tol=1e-10))
  probes <- list(truth=c(.10,.03,.02),low=c(.05,.015,.01),high=c(.3,.2,.05))
  for(probe in names(probes)) {
    p <- probes[[probe]]
    gp <- c(sA=p[1],sB=p[1],sAB=p[1],xA=p[3],xB=p[3],dA=p[2],dB=p[2])
    cp <- pars.geosse2classe(gp)
    tp <- c(t.1=z$recent.time,setNames(cp,paste0(names(cp),'.1')),setNames(cp,paste0(names(cp),'.2')))
    for(survival in c(FALSE,TRUE)) {
      ll.g <- g(gp,root=diversitree::ROOT.OBS,condition.surv=survival)
      ll.c <- cl.likelihood(tp,root=diversitree::ROOT.OBS,condition.surv=survival)
      audit[[length(audit)+1]] <- data.frame(id=z$id,probe=probe,survival=survival,difference=ll.c-ll.g)
    }
  }
}
control <- do.call(rbind,control)
audit <- do.call(rbind,audit)
write.csv(control,file.path(out,'geosse-estimate-comparison.csv'),row.names=FALSE)
write.csv(audit,file.path(out,'geosse-likelihood-audit.csv'),row.names=FALSE)

# Independently verify linked fits against unconstrained full parameter vectors. ----
link.audit <- list()
for(z in datasets) {
  raw <- make.classe.td(z$tree,z$tree$tip.state,k=3,n.epoch=2,strict=FALSE,control=list(tol=1e-10))
  best <- best.fits[[paste(z$id,'linked_classe')]]
  if(is.null(best)) next
  full <- z$truth
  for(rate in c('s','d','x')) {
    full[paste0(rate,'_recent')] <- best$par[paste0(rate,'_recent')]
    source <- if(z$links[rate]) paste0(rate,'_recent') else paste0(rate,'_old')
    full[paste0(rate,'_old')] <- best$par[source]
  }
  for(probe in c('truth','fit')) {
    p <- if(probe=='truth') z$truth else full
    values <- c(t.1=z$recent.time)
    epoch.index <- 0
    for(epoch in c('recent','old')) {
      epoch.index <- epoch.index+1
      g <- unname(p[paste0(c('s','s','s','x','x','d','d'),'_',epoch)])
      names(g) <- diversitree:::default.argnames.geosse()
      cp <- pars.geosse2classe(g)
      values <- c(values,setNames(cp,paste0(names(cp),'.',epoch.index)))
    }
    ll <- raw(values,root=diversitree::ROOT.OBS,condition.surv=TRUE)
    expected <- if(probe=='truth') best$truth.lnLik else best$lnLik
    link.audit[[length(link.audit)+1]] <- data.frame(id=z$id,probe=probe,difference=ll-expected)
  }
}
link.audit <- do.call(rbind,link.audit)
write.csv(link.audit,file.path(out,'parameter-link-audit.csv'),row.names=FALSE)
stopifnot(max(abs(link.audit$difference))<1e-4)

# Plots ----
p <- ggplot(estimates,aes(label,ratio,color=scenario))+
  geom_hline(yintercept=1,linetype=2)+geom_point(position=position_jitter(width=.12,height=0,seed=1))+
  facet_wrap(~scenario,scales='free_x',ncol=2)+scale_y_log10()+theme_classic()+
  theme(legend.position='none',axis.text.x=element_text(angle=40,hjust=1))+
  labs(x=NULL,y='Estimate / generating value',caption='Independent parameters only. Best converged fit; failed datasets are not plotted.')
ggsave(file.path(out,'recovery.png'),p,width=10,height=8,dpi=160)
ggsave(file.path(out,'recovery.pdf'),p,width=10,height=8)
p <- ggplot(estimates,aes(N,absolute.log10.error,color=epoch))+
  geom_point()+scale_x_log10()+facet_grid(rate~scenario)+theme_classic()+
  labs(x='Extant tips (log scale)',y='Absolute log10 estimation error',caption='Descriptive pilot: five trees per scenario; unequal sizes and different generating models.')
ggsave(file.path(out,'information-recovery.png'),p,width=12,height=6,dpi=160)
ggsave(file.path(out,'information-recovery.pdf'),p,width=12,height=6)
old <- estimates[estimates$epoch=='old',]
p <- ggplot(old,aes(boundary.lineages,absolute.log10.error,color=rate))+
  geom_point()+facet_wrap(~scenario)+theme_classic()+
  labs(x='Reconstructed lineages crossing epoch boundary',y='Absolute log10 estimation error')
ggsave(file.path(out,'old-information.png'),p,width=10,height=5,dpi=160)
ggsave(file.path(out,'old-information.pdf'),p,width=10,height=5)
p <- ggplot(control,aes(geosse,classe))+geom_abline(slope=1,intercept=0,linetype=2)+
  geom_point()+facet_wrap(~parameter,scales='free')+theme_classic()+
  labs(x='Classical GeoSSE estimate',y='Linked-epoch ClaSSE estimate')
ggsave(file.path(out,'geosse-comparison.png'),p,width=9,height=3,dpi=160)
ggsave(file.path(out,'geosse-comparison.pdf'),p,width=9,height=3)

# Report ----
summary <- list()
diagnostics <- list()
for(z in datasets) diagnostics[[z$id]] <- cbind(id=z$id,scenario=z$scenario,z$diagnostics)
diagnostics <- do.call(rbind,diagnostics)
write.csv(diagnostics,file.path(out,'tree-information.csv'),row.names=FALSE)
for(s in design$scenario_id) {
  cc <- checks[checks$scenario==s & checks$model=='linked_classe',]
  summary[[s]] <- data.frame(scenario=s,parameters=design$n_parameters[design$scenario_id==s],
    median.N=median(cc$N),median.crossing=median(cc$boundary.lineages),
    converged.starts=sum(cc$converged.starts),datasets.with.fit=sum(cc$converged.starts>0),
    all.three.agree=sum(cc$all.agree))
}
summary <- do.call(rbind,summary)
write.csv(summary,file.path(out,'scenario-summary.csv'),row.names=FALSE)
lines <- c('# Linked-rate recovery pilot','',
  'Five independent trees per scenario; three optimization starts per fit. All scenarios use 30 old + 6 recent time units. The boundary is fixed at age 6; this pilot does not vary boundary position.',
  'Simulation uses unchanged make.tree.classe.td(), uniformly sampled initial states A/B/AB, full sampling, and a 4000-taxon safety cap. Trees that stop early or lack a crown spanning both epochs are excluded; all attempts are logged.',
  'Inference imposes the scenario-specific links: 6/4/4/4/3 free parameters. Shared rates are listed once, not duplicated. ROOT.OBS, survival conditioning, ODE tolerance 1e-10, log-rate nlminb bounds 1e-10 to 100, three independently perturbed starts, 180-second limit per start.',
  'Convergence requires optimizer code zero and a finite likelihood. All-start agreement requires three converged starts within 1e-4 log-likelihood. It is not proof of identifiable parameters.',
  '', '## Design','', 's = speciation; d = range expansion; x = local extinction/range contraction.',
  '', '~~~',capture.output(print(design,row.names=FALSE)),'~~~',
  '', '## Convergence','',paste('Converged starts:',sum(starts$converged),'/',nrow(starts)),
  '', '~~~',capture.output(print(summary,row.names=FALSE)),'~~~',
  '', 'Nonconverged starts:','~~~',
  capture.output(print(starts[!starts$converged,c('id','model','start','message')],row.names=FALSE)),'~~~',
  '', '## Recovery','', '~~~',capture.output(print(recovery,row.names=FALSE)),'~~~',
  '', '## Classical GeoSSE comparison','',
  paste('Maximum matched-parameter log-likelihood difference:',max(abs(audit$difference))),
  paste('Maximum fitted log-likelihood difference:',max(abs(control$lnLik.difference))),
  paste('Maximum fitted parameter difference:',max(abs(control$classe-control$geosse))),
  paste('Maximum linked-versus-full ClaSSE likelihood difference:',max(abs(link.audit$difference))),
  paste('Linked fits below generating-parameter likelihood by >1e-4:',sum(estimates$lnLik<estimates$truth.lnLik-1e-4)),
  '', '## Information and limitations','',
  'tree-information.csv records extant tip counts by state, crown age, boundary age, crossing lineages, branch exposure and reconstructed branching events by epoch. It separately records full-history simulated speciation, terminal extinction, dispersal and contraction events; simulation tables and transition histories are retained in datasets.rds.',
  'The reconstructed root is counted as a branching event. A node exactly on the boundary belongs to the recent epoch.',
  'Five replicates per scenario cannot establish bias or coverage reliably. Sizes are outcomes, not controlled treatments. Plots of error against information are descriptive, not causal comparisons. No confidence intervals or coverage are estimated in this pilot.',
  'Models match the generating links; this is parameter recovery, not a test of model selection. No stochastic-mapping benchmark is performed here.',
  '', '## Figures','',
  '- [Recovery](recovery.pdf)',
  '- [Recovery versus tip number](information-recovery.pdf)',
  '- [Old-rate recovery versus crossing lineages](old-information.pdf)',
  '- [GeoSSE comparison](geosse-comparison.pdf)',
  '', '## Reproduce','',
  'Run R/12-linked-rate-pilot.R followed by R/13-report-linked-rate-pilot.R from the DivExtra repository root. Use a fresh output directory if changing settings; existing fits are reused.',
  '', '~~~',capture.output(sessionInfo()),'~~~')
writeLines(lines,file.path(out,'report.md'))
print(summary,row.names=FALSE)
print(recovery,row.names=FALSE)

