# Run after 06-fixed-time-pilot.R. Summarize every attempted inference.
library(ggplot2)
library(divextra)
out <- 'inst/benchmarks/geosse-parameter-recovery/results/fixed-time-pilot'
checks <- read.csv(file.path(out,'convergence.csv'))
starts <- read.csv(file.path(out,'starts.csv'))
est <- read.csv(file.path(out,'estimates.csv'))
design <- read.csv(file.path(out,'design.csv'))
attempts <- read.csv(file.path(out,'simulation-attempts.csv'))
fits <- readRDS(file.path(out,'fits.rds'))
datasets <- readRDS(file.path(out,'datasets.rds'))
est$error <- abs(log10(est$estimate/est$truth))

# Check that the selected fit improves on the likelihood at the true rates. ----
checks$truth.lnLik <- NA_real_
checks$best.minus.truth <- NA_real_
for (i in seq_along(datasets)) {
  z <- datasets[[i]]
  lik <- make.classe.td(z$tree,z$tree$tip.state,k=3,n.epoch=2,
    sampling.f=c(1,1,1),strict=FALSE,control=list(backend='gslode'))
  recent <- z$pars$classe[2,]
  old <- z$pars$classe[1,]
  true.par <- c(t.1=z$recent.time,
    setNames(recent,paste0(names(recent),'.1')),
    setNames(old,paste0(names(old),'.2')))
  checks$truth.lnLik[i] <- lik(true.par,root=diversitree::ROOT.OBS,condition.surv=TRUE)
  keep <- starts$id==z$id & starts$converged
  if (any(keep)) checks$best.minus.truth[i] <- max(starts$lnLik[keep])-checks$truth.lnLik[i]
}
write.csv(checks,file.path(out,'convergence.csv'),row.names=FALSE)

summary <- list()
for (s in design$scenario_id) {
  d <- checks[checks$scenario==s,]
  summary[[s]] <- data.frame(scenario=s,datasets=nrow(d),
    converged.starts=sum(d$converged.starts),all.three.agree=sum(d$all.starts.agree),
    median.N=median(d$N),median.cross=median(d$boundary.lineages),
    median.old.length=median(d$old.length))
}
summary <- do.call(rbind,summary)
write.csv(summary,file.path(out,'scenario-summary.csv'),row.names=FALSE)

recovery <- list()
for (s in design$scenario_id) {
  for (p in unique(est$parameter)) {
    d <- est[est$scenario==s & est$parameter==p,]
    recovery[[length(recovery)+1]] <- data.frame(scenario=s,parameter=p,n=nrow(d),
      truth=d$truth[1],median=median(d$estimate),mean=mean(d$estimate),
      minimum=min(d$estimate),maximum=max(d$estimate),
      median.absolute.log10.error=median(d$error),
      within.twofold=sum(d$estimate/d$truth>=.5 & d$estimate/d$truth<=2))
  }
}
recovery <- do.call(rbind,recovery)
write.csv(recovery,file.path(out,'recovery-summary.csv'),row.names=FALSE)

starts$delta.lnLik <- NA_real_
for (id in unique(starts$id)) {
  keep <- starts$id==id
  starts$delta.lnLik[keep] <- max(starts$lnLik[keep],na.rm=TRUE)-starts$lnLik[keep]
}
starts$id <- factor(starts$id,levels=checks$id)
p <- ggplot(starts,aes(id,factor(start),fill=log10(1+delta.lnLik)))+
  geom_tile(color='white')+scale_fill_viridis_c(na.value='grey80')+
  theme_classic()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  labs(x='Dataset',y='Starting point',fill='log10(1 + loss)',
    caption='Likelihood loss relative to the best returned start. Grey = missing/failed result.')
ggsave(file.path(out,'convergence.png'),p,width=10,height=3.5,dpi=180)
ggsave(file.path(out,'convergence.pdf'),p,width=10,height=3.5)

# Show the all-rates-change scenario on the original rate scale. ----
all.change <- est[est$scenario=='P05',]
all.change$epoch <- ifelse(grepl('_old',all.change$parameter),'Old','Recent')
all.change$rate <- sub('_.*','',all.change$parameter)
all.change$rate <- factor(all.change$rate,levels=c('s','d','x'),
  labels=c('Speciation','Dispersal','Contraction / extinction'))
truth <- unique(all.change[,c('epoch','rate','truth')])
p3 <- ggplot(all.change,aes(epoch,estimate,color=epoch))+
  geom_point(position=position_jitter(width=.1,height=0,seed=1),size=2.5)+
  geom_point(data=truth,aes(epoch,truth),inherit.aes=FALSE,shape=3,size=5)+
  facet_wrap(~rate,scales='free_y')+theme_classic()+
  scale_color_manual(values=c(Old='#0072B2',Recent='#D55E00'))+
  theme(legend.position='none')+labs(x=NULL,y='Estimated rate',
    caption='P05: all three rates change. Five datasets; black crosses show the generating rates.')
ggsave(file.path(out,'all-rates-change.png'),p3,width=9,height=4,dpi=180)
ggsave(file.path(out,'all-rates-change.pdf'),p3,width=9,height=4)

lines <- c('# Fixed-time GeoSSE recovery pilot', '',
  'This pilot uses five independently simulated datasets in each of five scenarios.',
  'Each dataset was fitted with the six-rate episodic model from three starting vectors.',
  'Starting vectors use eps = 0.1, 0.4, 0.8 and independent multiplicative perturbations between 0.5 and 2.',
  'Inference uses ROOT.OBS and survival conditioning. Each optimizer has a 90-second wall-time limit; three processes run concurrently.',
  'Optimizer code zero and finite likelihood are required to select a best converged estimate.',
  'Likelihood agreement means all three starts converge and their log-likelihood range is below 0.0001.',
  'Parameter relative spread is also recorded: equal likelihoods do not establish equal parameters or identifiability.',
  '', '## Duration calibration', '',
  'For each scenario, 20 independent attempts were simulated for each of ten duration combinations.',
  'Candidates required at least eight usable trees and no taxon-limit truncations in calibration.',
  'The score targeted 30 reconstructed lineages at the boundary and approximately 200 extant tips.',
  'The chosen durations were frozen before independent recovery simulations. No filtering on realized tip count or boundary lineage count was applied.',
  'Trees lacking a reconstructed crown spanning both epochs were excluded; this is additional design conditioning.',
  'A 3000-taxon cap protects computation. This small calibration cannot guarantee caps will never occur.',
  '', '```',capture.output(print(design[,c('scenario_id','old.time','recent.time')],row.names=FALSE)),
  '```', '', '## Completed fits', '',
  paste('Converged starts:',sum(starts$converged),'of',nrow(starts)),
  paste('Datasets with three agreeing converged starts:',sum(checks$all.starts.agree),'of',nrow(checks)),
  paste('Datasets whose best converged likelihood is below truth by more than 0.0001:',
    sum(checks$best.minus.truth < -1e-4,na.rm=TRUE)),
  '', '```',capture.output(print(summary,row.names=FALSE)), '```',
  '', '## All three rates change: P05', '',
  'Old rates: s=0.10, d=0.03, x=0.02. Recent rates: s=0.30, d=0.20, x=0.15.',
  '', '```',capture.output(print(recovery[recovery$scenario=='P05',],row.names=FALSE)), '```',
  '', '## Simulation exclusions', '', '```',
  capture.output(print(table(attempts$scenario,attempts$status))), '```',
  '', '## Limits', '',
  'Five trees per scenario support a diagnostic pilot, not precise bias or coverage estimates.',
  'Recovery plots show individual estimates. No 95% confidence interval is claimed from these five replicates.',
  'The scenario durations differ, so comparisons with the earlier taxon-stopped pilot are descriptive rather than paired experiments.',
  'Convergence codes, agreement between starts, parameter recovery, and statistical identifiability are distinct checks.',
  'Simulation cap exclusions and optimizer timeouts must be included in production reports, not silently discarded.',
  '', '## Pilot interpretation and next steps', '',
  'Recent speciation and dispersal generally recovered better than old dispersal and contraction/extinction.',
  'P05 old dispersal and old contraction estimates were outside a factor of two of truth in all five datasets.',
  'Increased old-epoch history did not by itself resolve recovery of these parameters.',
  'Likelihood improvement over truth rules out simply failing to reach the generating likelihood; it does not prove the global optimum or validate the entire likelihood implementation.',
  'Parameter spread can be large for rates near zero despite almost identical likelihoods; see the saved fit vectors.',
  'Before a large run, profile old dispersal and contraction for representative P05 datasets and check likelihood trade-offs and numerical accuracy.',
  'Follow with larger independent replicates spanning old-epoch information levels; report all convergence failures and boundary estimates.',
  '', '## Reproduction', '',
  'Run R/06-fixed-time-pilot.R, followed by R/07-report-fixed-time-pilot.R, from the DivExtra repository root.',
  'datasets.rds contains the simulated trees and histories; fits.rds includes the initial vectors and full fit objects.',
  '', '```',capture.output(sessionInfo()),'```')
writeLines(lines,file.path(out,'report.md'))
print(summary,row.names=FALSE)
print(recovery[recovery$scenario=='P05',],row.names=FALSE)
