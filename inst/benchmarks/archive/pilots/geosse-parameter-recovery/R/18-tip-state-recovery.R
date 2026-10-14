# Existing epoch-grid pilot only: no simulations, filtering or refitting. ----
library(ggplot2)
base <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/epoch-grid-pilot'
out <- file.path(base,'tip-state-diagnostic')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tips <- read.csv(file.path(base,'tree-information.csv'))
estimates <- read.csv(file.path(base,'regression-data.csv'))
stopifnot(tips$n.A+tips$n.B+tips$n.AB==tips$N)
tips$min.count <- pmin(tips$n.A,tips$n.B,tips$n.AB)
tips$min.fraction <- tips$min.count/tips$N
tips$AB.fraction <- tips$n.AB/tips$N
tips$all.present <- tips$min.count>=1
tips$five.each <- tips$min.count>=5
write.csv(tips,file.path(out,'tip-state-counts.csv'),row.names=FALSE)
counts <- aggregate(cbind(n.A,n.B,n.AB,min.count)~scenario,tips,median)
counts$all.present <- as.integer(tapply(tips$all.present,tips$scenario,sum)[counts$scenario])
counts$five.each <- as.integer(tapply(tips$five.each,tips$scenario,sum)[counts$scenario])
write.csv(counts,file.path(out,'state-summary.csv'),row.names=FALSE)
d <- merge(estimates,tips[,c('id','min.count','min.fraction','AB.fraction','all.present','five.each')],by='id')
d$within.twofold <- d$ratio>=.5 & d$ratio<=2
d$near.zero <- d$estimate<1e-5

# Threshold summaries are descriptive, not independent treatment groups. ----
summary <- list()
for(label in unique(d$label)) for(rule in c('all','five.each')) {
  keep <- d$label==label
  if(rule=='five.each') keep <- keep & d$five.each
  z <- d[keep,]
  summary[[length(summary)+1]] <- data.frame(parameter=label,rule=rule,n=nrow(z),
    within.twofold=sum(z$within.twofold),near.zero=sum(z$near.zero),
    median.bounded.error=median(z$bounded.error))
}
summary <- do.call(rbind,summary)
write.csv(summary,file.path(out,'threshold-recovery.csv'),row.names=FALSE)

# Continuous associations adjusted for tree size and scenario-boundary cell. ----
regression <- list()
for(label in unique(d$label)) for(predictor in c('min.fraction','AB.fraction')) {
  z <- d[d$label==label,]
  z$representation <- z[[predictor]]/.10
  for(outcome in c('absolute.log10.error','bounded.error')) {
    z$error <- z[[outcome]]
    fit <- lm(error~representation+log2(N)+factor(scenario),data=z)
    X <- model.matrix(fit)
    if(fit$rank<ncol(X) || df.residual(fit)<5) next
    bread <- solve(crossprod(X))
    weighted <- X*as.numeric(residuals(fit)/(1-hatvalues(fit)))
    variance <- bread %*% crossprod(weighted) %*% bread
    se <- sqrt(variance['representation','representation'])
    slope <- coef(fit)['representation']
    half <- qt(.975,df.residual(fit))*se
    regression[[length(regression)+1]] <- data.frame(parameter=label,predictor=predictor,
      outcome=outcome,n=nrow(z),slope.per.ten.percentage.points=unname(slope),
      lower=unname(slope-half),upper=unname(slope+half))
  }
}
regression <- do.call(rbind,regression)
write.csv(regression,file.path(out,'state-recovery-regressions.csv'),row.names=FALSE)
old <- d[d$label=='x_old',]
p <- ggplot(old,aes(min.count,estimate,color=scenario))+
  geom_hline(yintercept=.02,linetype=2)+geom_point(size=2)+
  geom_vline(xintercept=5,linetype=3)+scale_x_log10()+scale_y_log10()+theme_classic()+
  labs(x='Smallest tip-state count (log scale)',y='Estimated old contraction (log scale)',
    caption='Each point is one existing tree. Dashed: true rate 0.02. Dotted: five tips per state.')
ggsave(file.path(out,'old-contraction-tip-states.png'),p,width=10,height=5,dpi=160)
ggsave(file.path(out,'old-contraction-tip-states.pdf'),p,width=10,height=5)
lines <- c('# Tip-state representation and recovery','',
  '## Main findings','',
  'All 75 trees already contain A, B and AB, and 72 meet the five-tips-per-state rule. Missing-state trees cannot explain the existing recovery failures.',
  'For old contraction, 21/30 estimates are below 1e-5. Applying the five-tip rule leaves 29 trees and all 21 near-zero estimates; recovery within a factor of two changes only from 5/30 to 5/29. Thus this rule does not resolve the problem in the existing sample.',
  'After controlling tree size and scenario-boundary cell, old-contraction association intervals include zero for both state-balance predictors and both error definitions. The sample does not establish absence of an effect, but supplies no clear evidence of improved old-contraction recovery with more balanced present-day states.',
  'Across all 36 exploratory regressions, two nominal intervals exclude zero, both in the opposite direction (more AB associated with worse log-scale error for recent/shared contraction). Neither is reproduced with bounded error. These are not multiplicity-adjusted findings and do not support choosing a state threshold.',
  '', '## Procedure','',
  'Existing 75-tree epoch-grid pilot only. No new simulations, rejections or inference fits. State counts were verified to sum to tip number.',
  paste('All three states present:',sum(tips$all.present),'of',nrow(tips)),
  paste('At least five tips per state:',sum(tips$five.each),'of',nrow(tips)),
  '', '## Median state counts by cell','',
  'These are separate column medians and need not sum to median tree size.',
  '~~~',capture.output(print(counts,row.names=FALSE)),'~~~',
  '', '## Trees failing five tips per state','',
  '~~~',capture.output(print(tips[!tips$five.each,c('id','N','n.A','n.B','n.AB')],row.names=FALSE)),'~~~',
  '', '## Recovery with and without the proposed five-tip rule','',
  'Summaries are by independently estimated parameter. The two groups overlap; they are not randomized treatments.',
  '~~~',capture.output(print(summary,row.names=FALSE)),'~~~',
  '', '## Continuous association checks','',
  'Separate models for each parameter and state predictor: error ~ state proportion + log2(tip number) + scenario-boundary cell. Slopes are per 10 percentage-point increase; HC3 95% intervals. Outcomes are absolute log10 estimate/truth error and bounded symmetric error. Negative slopes mean better recovery.',
  'No missing-state trees exist here, and only three trees fail the five-tip rule. This dataset cannot estimate the benefit of requiring all states or give a reliable threshold-treatment comparison. Continuous results are exploratory, with no multiple-testing adjustment. Five trees per cell and correlated predictors limit precision.',
  '~~~',capture.output(print(regression,row.names=FALSE)),'~~~',
  '', '[Old contraction versus smallest state count](old-contraction-tip-states.pdf)')
writeLines(lines,file.path(out,'report.md'))
print(summary,row.names=FALSE)
print(regression[regression$parameter=='x_old',],row.names=FALSE)
