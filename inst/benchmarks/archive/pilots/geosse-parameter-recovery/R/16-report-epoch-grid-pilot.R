# Report recovery, convergence and the classical GeoSSE equivalence check.
library(divextra)
library(ggplot2)
out <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/epoch-grid-pilot'
datasets <- readRDS(file.path(out,'datasets.rds'))
fits <- readRDS(file.path(out,'fits.rds'))
design <- read.csv(file.path(out,'design.csv'))
stopifnot(length(datasets)==75, length(fits)==270)
for(f in fits) stopifnot(!is.null(f))
for(z in datasets) stopifnot(z$diagnostics$N>=50,z$diagnostics$N<=1000,
  z$diagnostics$old.nodes>=10,z$diagnostics$recent.nodes>=10,
  z$diagnostics$boundary.lineages>=10)
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
estimates$rate.scenario <- design$rate.scenario[match(estimates$scenario,design$scenario_id)]
estimates$fraction <- design$fraction[match(estimates$scenario,design$scenario_id)]
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
for(z in datasets) if(z$rate.scenario=='constant') {
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
  facet_wrap(~scenario,scales='free_x',ncol=3)+scale_y_log10()+theme_classic()+
  theme(legend.position='none',axis.text.x=element_text(angle=40,hjust=1))+
  labs(x=NULL,y='Estimate / generating value',caption='Independent parameters only. Best converged fit; failed datasets are not plotted.')
ggsave(file.path(out,'recovery.png'),p,width=13,height=15,dpi=160)
ggsave(file.path(out,'recovery.pdf'),p,width=13,height=15)
p <- ggplot(estimates,aes(N,absolute.log10.error,color=epoch))+
  geom_point()+scale_x_log10()+facet_grid(rate~rate.scenario)+theme_classic()+
  labs(x='Extant tips (log scale)',y='Absolute log10 estimation error',caption='Five trees per scenario-boundary cell; unequal sizes and calibrated durations. Descriptive only.')
ggsave(file.path(out,'information-recovery.png'),p,width=12,height=6,dpi=160)
ggsave(file.path(out,'information-recovery.pdf'),p,width=12,height=6)
old <- estimates[estimates$epoch=='old',]
p <- ggplot(old,aes(boundary.lineages,absolute.log10.error,color=rate))+
  geom_point()+facet_wrap(~scenario)+theme_classic()+
  labs(x='Reconstructed lineages crossing epoch boundary',y='Absolute log10 estimation error')
ggsave(file.path(out,'old-information.png'),p,width=13,height=10,dpi=160)
ggsave(file.path(out,'old-information.pdf'),p,width=13,height=10)
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
tree.plot <- rbind(data.frame(id=diagnostics$id,scenario=diagnostics$scenario,
  metric='Extant tips',value=diagnostics$N),
  data.frame(id=diagnostics$id,scenario=diagnostics$scenario,
  metric='Crossing lineages',value=diagnostics$boundary.lineages))
tree.plot$rate.scenario <- design$rate.scenario[match(tree.plot$scenario,design$scenario_id)]
tree.plot$boundary <- factor(design$fraction[match(tree.plot$scenario,design$scenario_id)],
  levels=c(.05,.25,.5),labels=c('5%','25%','50%'))
p <- ggplot(tree.plot,aes(boundary,value))+
  geom_point(position=position_jitter(width=.1,height=0,seed=1))+
  facet_grid(metric~rate.scenario,scales='free_y')+theme_classic()+
  labs(x='Boundary fraction of calibrated simulation duration',y=NULL,
    caption='Five accepted trees in every cell; durations calibrated separately. All tips between 50 and 1000.')
ggsave(file.path(out,'tree-information.png'),p,width=12,height=5,dpi=160)
ggsave(file.path(out,'tree-information.pdf'),p,width=12,height=5)
for(s in design$scenario_id) {
  cc <- checks[checks$scenario==s & checks$model=='linked_classe',]
  summary[[s]] <- data.frame(scenario=s,parameters=design$n_parameters[design$scenario_id==s],
    median.N=median(cc$N),median.crossing=median(cc$boundary.lineages),
    converged.starts=sum(cc$converged.starts),datasets.with.fit=sum(cc$converged.starts>0),
    all.three.agree=sum(cc$all.agree))
}
summary <- do.call(rbind,summary)
attempts <- read.csv(file.path(out,'simulation-attempts.csv'))
summary$attempts <- as.integer(table(factor(attempts$scenario,levels=summary$scenario)))
summary$acceptance <- 5/summary$attempts
summary$total.duration <- design$total.time[match(summary$scenario,design$scenario_id)]
write.csv(summary,file.path(out,'scenario-summary.csv'),row.names=FALSE)
# Associations with information, controlling each scenario-boundary cell. ----
# Each model has one independently estimated parameter per tree.
# No pooling of multiple parameters as independent observations.
regression.data <- merge(estimates,diagnostics[,c('id','recent.length')],by='id')
regression.data$branch.exposure <- regression.data$old.length+regression.data$recent.length
regression.data$branch.exposure[regression.data$epoch=='old'] <- regression.data$old.length[regression.data$epoch=='old']
regression.data$branch.exposure[regression.data$epoch=='recent'] <- regression.data$recent.length[regression.data$epoch=='recent']
regression.data$bounded.error <- abs(regression.data$estimate-regression.data$truth)/(regression.data$estimate+regression.data$truth)
regression.rows <- list()
for(label in unique(regression.data$label)) {
  d <- regression.data[regression.data$label==label,]
  for(predictor in c('N','boundary.lineages','branch.exposure')) {
    d$information <- log2(d[[predictor]])
    for(outcome in c('absolute.log10.error','bounded.error')) {
      d$error <- d[[outcome]]
      fit <- lm(error ~ information + factor(scenario),data=d)
      X <- model.matrix(fit)
      if(fit$rank<ncol(X) || df.residual(fit)<5) next
      # HC3 heteroscedasticity-robust standard errors.
      bread <- solve(crossprod(X))
      weighted <- X*as.numeric(residuals(fit)/(1-hatvalues(fit)))
      variance <- bread %*% crossprod(weighted) %*% bread
      se <- sqrt(variance['information','information'])
      slope <- coef(fit)['information']
      half <- qt(.975,df.residual(fit))*se
      regression.rows[[length(regression.rows)+1]] <- data.frame(parameter=label,
        predictor=predictor,outcome=outcome,n=nrow(d),cells=length(unique(d$scenario)),
        slope.per.doubling=unname(slope),lower=unname(slope-half),upper=unname(slope+half))
    }
  }
}
regression <- do.call(rbind,regression.rows)
write.csv(regression,file.path(out,'information-regressions.csv'),row.names=FALSE)
write.csv(regression.data,file.path(out,'regression-data.csv'),row.names=FALSE)
p <- ggplot(regression,aes(parameter,slope.per.doubling,ymin=lower,ymax=upper))+
  geom_hline(yintercept=0,linetype=2)+geom_pointrange()+
  facet_grid(outcome~predictor,scales='free_y')+theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust=1))+
  labs(x=NULL,y='Change in error per doubling (HC3 95% interval)',
    caption='Cell-adjusted exploratory associations; separate predictors, no causal interpretation or multiplicity correction.')
ggsave(file.path(out,'information-regressions.png'),p,width=12,height=7,dpi=160)
ggsave(file.path(out,'information-regressions.pdf'),p,width=12,height=7)

lines <- c('# Calibrated epoch-grid recovery pilot','',
  '## Main findings','',
  'All 15 scenario-boundary combinations are represented by five qualifying trees. Of 270 optimization starts, 269 converged and one timed out; every tree has a converged fit. All three ClaSSE starts agree in likelihood for 71 of 75 trees. No returned estimate hits the upper parameter bound.',
  'Recent speciation and dispersal medians are generally close to truth. Old contraction remains difficult, often approaching zero. In the all-rates-change scenario its median is near zero at all three boundaries. The information rules therefore do not guarantee recovery.',
  'The cell-adjusted regressions suggest smaller recent-speciation errors with more tips, crossing lineages and recent branch exposure under both error metrics. Shared speciation improves with tip count, and shared dispersal with crossing lineages/branch exposure. Contraction slopes remain uncertain. These are exploratory patterns across 54 related fits, not multiplicity-adjusted discoveries.',
  'At 50%, old speciation medians are elevated in the all-rates-change and speciation-only scenarios (about 0.169 and 0.168 versus truth 0.10). Qualifying trees are rare in these cells, so results must be interpreted conditional on selection. The pilot cannot separate selection effects from other sources of error.',
  '', '## Procedure','',
  'Five independent trees per scenario-boundary cell, 15 cells, 75 trees. Boundaries are 5%, 25%, 50% of simulation duration measured backwards from present, NOT reconstructed crown age. Total durations are calibrated per cell using 30 independent attempts at each of 18,24,30,36,42,48 units; highest acceptance wins, ties favor median size nearest 200.',
  'Simulation uses unchanged make.tree.classe.td(), uniformly sampled initial states A/B/AB, full sampling, and a 1000-taxon safety cap. Final trees require 50--1000 tips, at least 10 crossing lineages and 10 reconstructed branching events per epoch. All attempts are logged. Calibration never uses recovery quality; final trees use fresh seeds.',
  'Inference imposes the scenario-specific links: 6/4/4/4/3 free parameters. Shared rates are listed once, not duplicated. ROOT.OBS, survival conditioning, ODE tolerance 1e-10, log-rate nlminb bounds 1e-10 to 100, three independently perturbed starts, 180-second limit per start.',
  'Convergence requires optimizer code zero and a finite likelihood. All-start agreement requires three converged starts within 1e-4 log-likelihood. It is not proof of identifiable parameters.',
  '', '## Design','', 's = speciation; d = range expansion; x = local extinction/range contraction.',
  '', '~~~',capture.output(print(design,row.names=FALSE)),'~~~',
  '', '## Convergence','',paste('Converged starts:',sum(starts$converged),'/',nrow(starts)),
  paste('Returned starts at upper parameter bound:',sum(starts$bound.hit,na.rm=TRUE)),
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
  'Five replicates per cell cannot establish bias or coverage reliably. Sizes are outcomes, not controlled treatments. Plots of error against information are descriptive, not causal comparisons. No parameter-recovery confidence intervals or coverage are estimated in this pilot; regression slope intervals are reported separately.',
  'The all-rates-change 50% cell had zero accepted trees in its 30-attempt calibration batches. The deterministic fallback selected the shortest candidate duration (18); the independent final sample required 500 attempts to obtain five qualifying trees. This very low acceptance must not be hidden when reporting recovery.',
  'Models match the generating links; this is parameter recovery, not a test of model selection. No stochastic-mapping benchmark is performed here.',
  '', '## Information regressions','',
  'Separate regressions for each independent parameter and each predictor: error ~ log2(information) + scenario-boundary cell fixed effects. Branch exposure is epoch-specific for epoch-specific rates and total exposure for shared rates. This controls duration and rate settings through cell effects. Crossing counts are not fitted jointly with tip counts or branch exposure because they are correlated.',
  'Outcomes: absolute log10 estimate/truth error, and bounded symmetric error abs(estimate-truth)/(estimate+truth). Negative slopes mean better recovery with more information. HC3 robust 95% intervals are exploratory, not adjusted for multiple comparisons. Sparse cells, selected trees, boundary estimates and failed fits limit inference; no causal conclusion.',
  '', '~~~',capture.output(print(regression,row.names=FALSE)),'~~~',
  '', '## Figures','',
  '- [Information regressions](information-regressions.pdf)',
  '- [Tree sizes and crossing lineages](tree-information.pdf)',
  '- [Recovery](recovery.pdf)',
  '- [Recovery versus tip number](information-recovery.pdf)',
  '- [Old-rate recovery versus crossing lineages](old-information.pdf)',
  '- [GeoSSE comparison](geosse-comparison.pdf)',
  '', '## Reproduce','',
  'Run R/15-epoch-grid-pilot.R followed by R/16-report-epoch-grid-pilot.R from the DivExtra repository root. Use a fresh output directory if changing settings; existing fits are reused.',
  '', '~~~',capture.output(sessionInfo()),'~~~')
writeLines(lines,file.path(out,'report.md'))
print(summary,row.names=FALSE)
print(recovery,row.names=FALSE)


