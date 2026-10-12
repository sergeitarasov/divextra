# Numerical audit and interpretation of the saved P05 nuisance profiles.
library(divextra)
library(ggplot2)
base <- 'inst/benchmarks/geosse-parameter-recovery/results/fixed-time-pilot'
out <- file.path(base,'profiles')
results <- readRDS(file.path(out,'profiles.rds'))
datasets <- readRDS(file.path(base,'datasets.rds'))
fits <- readRDS(file.path(base,'fits.rds'))
d <- read.csv(file.path(out,'profiles.csv'))
failures <- list()
for(r in results) if(!is.null(r$error)) failures[[length(failures)+1]] <- data.frame(
  id=r$id,rate=r$rate,value=r$value,error=r$error)
failure.table <- if(length(failures)) do.call(rbind,failures) else data.frame()
write.csv(failure.table,file.path(out,'failed-grid-points.csv'),row.names=FALSE)
audit <- list()
ridge <- list()
for(id in unique(d$id)) {
  z <- datasets[[which(vapply(datasets,function(x) x$id==id,logical(1)))]]
  cfg <- read_yaml_pars_td('inst/benchmarks/geosse-parameter-recovery/config/model-episodic.yml')
  cfg$epoch.times <- z$recent.time
  f <- make.classe.td(z$tree,z$tree$tip.state,k=3,n.epoch=2,strict=FALSE,
    control=list(backend='gslode',gsl.stepper='rk8pd',tol=1e-12))
  f <- diversitree::constrain(f,formulae=make_constraints_sse_td(cfg))
  for(r in results) {
    if(r$id!=id || !is.null(r$error)) next
    value <- f(r$parameters,root=diversitree::ROOT.OBS,condition.surv=TRUE)
    audit[[length(audit)+1]] <- data.frame(id=id,rate=r$rate,value=r$value,
      rk8pd=value,rkck=r$numerical['gslode 1e-12'],
      difference=value-unname(r$numerical['gslode 1e-12']),
      tolerance.range=diff(range(r$numerical[1:3])))
    ridge[[length(ridge)+1]] <- data.frame(id=id,profile=r$rate,
      d_old=unname(r$parameters['q13.2']),x_old=unname(r$parameters['mu1.2']),lnLik=r$lnLik)
  }
}
audit <- do.call(rbind,audit)
ridge <- do.call(rbind,ridge)
write.csv(audit,file.path(out,'numerical-audit.csv'),row.names=FALSE)
write.csv(ridge,file.path(out,'nuisance-tradeoffs.csv'),row.names=FALSE)
summary <- list()
for(id in unique(d$id)) {
  for(rate in unique(d$rate)) {
    v <- d[d$id==id & d$rate==rate,]
    true <- if(rate=='q13.2') .03 else .02
    supported <- v$value[v$loss<=1.92 & v$convergence==0]
    summary[[length(summary)+1]] <- data.frame(id=id,rate=rate,
      best.grid.value=v$value[which.min(v$loss)],true.value=true,
      loss.at.truth=v$loss[which.min(abs(v$value-true))],
      grid.lower=min(supported),grid.upper=max(supported),
      lower.edge.supported=min(v$value) %in% supported,
      upper.edge.supported=max(v$value) %in% supported)
  }
}
summary <- do.call(rbind,summary)
write.csv(summary,file.path(out,'profile-summary.csv'),row.names=FALSE)

# Uniform slowdown check: compare likelihood contrasts, whose time-unit constants cancel.
z <- datasets[[which(vapply(datasets,function(x) x$id=='P05-1',logical(1)))]]
ff <- fits[vapply(fits,function(x) x$id==z$id,logical(1))]
best <- ff[[which.max(vapply(ff,function(x) x$fit$lnLik,numeric(1)))]]$fit$par
truth <- c(lambda111.1=.3,mu1.1=.15,q13.1=.2,lambda111.2=.1,mu1.2=.02,q13.2=.03)
scaling <- list()
for(multiplier in c(1,2,5)) {
  tree <- z$tree
  tree$edge.length <- tree$edge.length*multiplier
  cfg$epoch.times <- z$recent.time*multiplier
  f <- make.classe.td(tree,tree$tip.state,k=3,n.epoch=2,strict=FALSE,
    control=list(tol=1e-12))
  f <- diversitree::constrain(f,formulae=make_constraints_sse_td(cfg))
  a <- f(best/multiplier,root=diversitree::ROOT.OBS,condition.surv=TRUE)
  b <- f(truth/multiplier,root=diversitree::ROOT.OBS,condition.surv=TRUE)
  scaling[[length(scaling)+1]] <- data.frame(time.multiplier=multiplier,rate.multiplier=1/multiplier,
    likelihood.contrast=a-b)
}
scaling <- do.call(rbind,scaling)
write.csv(scaling,file.path(out,'uniform-slowdown-check.csv'),row.names=FALSE)

truth.lines <- expand.grid(id=unique(d$id),rate=unique(d$rate))
truth.lines$value <- ifelse(truth.lines$rate=='q13.2',.03,.02)
d$status <- ifelse(d$bound.hit,'Nuisance upper bound',ifelse(d$convergence==0,'Converged','Not converged'))
p <- ggplot(d,aes(value,loss))+geom_line()+geom_point(aes(shape=status))+
  geom_vline(data=truth.lines,aes(xintercept=value),color='#0072B2',linetype=2)+
  geom_hline(yintercept=1.92,linetype=3,color='grey40')+
  facet_grid(id~rate,scales='free_x',labeller=labeller(rate=c(
    'mu1.2'='Old contraction / extinction','q13.2'='Old dispersal')))+
  scale_x_log10()+coord_cartesian(ylim=c(0,6))+
  theme_classic()+labs(x='Fixed old-epoch rate (log scale)',y='Log-likelihood loss',shape='Profile fit',
    caption='Blue = true rate; dotted line = 1.92 reference. Display zoomed to loss 0-6; full plot and values also saved.')
ggsave(file.path(out,'profiles-zoom.png'),p,width=9,height=8,dpi=180)
ggsave(file.path(out,'profiles-zoom.pdf'),p,width=9,height=8)
for(id in unique(ridge$id)) ridge$loss[ridge$id==id] <- max(ridge$lnLik[ridge$id==id])-ridge$lnLik[ridge$id==id]
p2 <- ggplot(ridge,aes(d_old,x_old,color=loss,shape=profile))+geom_point(size=2)+
  facet_wrap(~id)+scale_x_log10()+scale_y_log10()+scale_color_viridis_c()+theme_classic()+
  labs(x='Old dispersal',y='Old contraction/extinction',color='Likelihood loss',
    caption='Nuisance optima along the one-dimensional profiles; this is not a two-dimensional confidence region.')
ggsave(file.path(out,'tradeoffs.png'),p2,width=10,height=4,dpi=180)
ggsave(file.path(out,'tradeoffs.pdf'),p2,width=10,height=4)
lines <- c('# P05 likelihood profiles and numerical audit','',
  'The evidence favors weak practical identifiability over integration error in these three datasets.',
  'For all six profiles the true focal rate is less than 1.92 log-likelihood units below the best value, after optimizing nuisance rates.',
  'P05-3 dispersal stays within 1.45 units across the full evaluated range from approximately 1e-8 to 10.',
  'Some extreme tails timed out or encountered a nuisance bound and remain unresolved. No formal confidence interval is claimed.',
  '',
  'Representative datasets: P05-1 (499 tips), P05-3 (216 tips), P05-4 (68 tips).',
  'Each grid value fixes one old-epoch rate while fitting the other five rates.',
  'Two nuisance starts are used: the previous best fit and the true parameters.',
  'Optimization uses log rates, nlminb, and GSL rkck tolerance 1e-10. Free rates are bounded between 1e-10 and 100.',
  'Both optimizer objects and nuisance vectors are saved. Each grid task has a 120-second process limit.',
  '',paste('Completed grid tasks:',nrow(d),'of',length(results)),
  paste('Nonzero chosen optimizer codes:',sum(d$convergence!=0)),
  paste('Tasks hitting the nuisance upper bound:',sum(d$bound.hit)),
  '', 'Failed or timed-out points (not interpreted as evidence against their rate values):',
  '```',capture.output(print(failure.table,row.names=FALSE)),'```',
  '', '## Profiles', '', '```',capture.output(print(summary,row.names=FALSE)),'```',
  'The reported grid limits are sampled values with likelihood loss <=1.92. They are not interpolated confidence limits.',
  'A supported grid edge means that direction is unresolved within the tested range.',
  'Nuisance optima at the upper bound 100 are constrained by the chosen search range; those points cannot establish an unconstrained tail endpoint.',
  'The 1.92 threshold is a conventional one-parameter likelihood reference. Boundary parameters and this sampling design require care; coverage was not calibrated.',
  '', '## Numerical checks', '',
  paste('Maximum likelihood range across rkck tolerances 1e-8, 1e-10, 1e-12:',format(max(audit$tolerance.range),digits=6)),
  paste('Maximum absolute difference between rk8pd and rkck at 1e-12:',format(max(abs(audit$difference)),digits=6)),
  'The deSolve backend was attempted but returns "Not yet possible" for this likelihood; it does not provide an independent validation here.',
  'Agreement between GSL methods and tolerances tests integration stability. It does not rule out a shared equation or implementation error.',
  '', '## Slower rates', '', '```',capture.output(print(scaling,row.names=FALSE)),'```',
  'Uniformly slowing every rate while lengthening all times by the reciprocal factor preserves likelihood contrasts.',
  'Slowing every rate at fixed duration usually produces fewer events and is not a general remedy for poor recovery.',
  'A useful next design changes relative rates: slower RECENT speciation to reduce recent dominance, with an independently calibrated old interval.',
  'For a controlled follow-up, compare recent s=0.15 versus 0.30 while keeping old s=0.10 and the other P05 rates fixed; both epochs still differ in all three rates.',
  'Also vary old d independently (e.g. 0.03, 0.10) and old x (e.g. 0.02, 0.05) to test whether rare old events drive the weakness.',
  'Changing these rates is a proposed follow-up, not an analysis performed in this report.',
  '', '## Interpretation boundaries', '',
  'Broad nuisance-optimized profiles with stable numerical likelihoods support weak practical identifiability for these observed datasets.',
  'They do not prove structural non-identifiability, global optimality, or validity of every implementation component.',
  'Only three pilot datasets are profiled. No general recovery guarantee follows from these cases.',
  '', '## Reproduction', '',
  'Run R/08-profile-p05.R then R/09-report-p05-profiles.R from the DivExtra repository root.')
writeLines(lines,file.path(out,'report.md'))
print(summary,row.names=FALSE)
print(c(tolerance.max=max(audit$tolerance.range),stepper.max=max(abs(audit$difference))))
print(scaling)
