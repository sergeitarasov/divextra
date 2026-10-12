# Simulation-only screening of the proposed information rules. ----
# Run from the DivExtra repository root. No parameter inference is performed.
library(divextra)
source('inst/benchmarks/geosse-parameter-recovery/R/benchmark-functions.R')
out <- 'inst/benchmarks/geosse-parameter-recovery/results/information-screen'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
design <- read.csv('inst/benchmarks/geosse-parameter-recovery/config/scenarios-linked.csv')
fractions <- c(.25, .50, .75)
n.attempts <- 50
rows <- list()

for(s in seq_len(nrow(design))) {
  parameters <- make_simulation_parameters(design[s,])
  for(f in seq_along(fractions)) {
    recent.time <- 36*fractions[f]
    old.time <- 36-recent.time
    for(attempt in seq_len(n.attempts)) {
      seed <- 14100000+s*10000+f*1000+attempt
      set.seed(seed)
      root <- sample(1:3, 1)
      row <- data.frame(scenario=design$scenario_id[s], boundary.fraction=fractions[f],
        old.duration=old.time, recent.duration=recent.time, attempt=attempt,
        seed=seed, root=root, status='complete', N=NA_integer_, crossing=NA_integer_,
        old.nodes=NA_integer_, recent.nodes=NA_integer_, crown.age=NA_real_,
        old.length=NA_real_, recent.length=NA_real_, accepted=FALSE)
      tb <- make.tree.classe.td(parameters$classe, k=3, max.taxa=4000,
        max.t1=old.time, max.t2=36, x0=root, single.lineage=TRUE)
      if(is.null(tb)) row$status <- 'no_tree'
      if(row$status=='complete' && !any(!tb$split & !tb$extinct)) row$status <- 'no_tree'
      if(row$status=='complete' && (abs(attr(tb,'t.epoch1')-old.time)>1e-8 ||
        abs(attr(tb,'t.total')-36)>1e-8)) row$status <- 'taxon_cap'
      if(row$status=='complete') {
        tree <- tryCatch(table2tree(tb), error=function(e) NULL)
        if(is.null(tree) || ape::Ntip(tree)<2) row$status <- 'no_tree'
      }
      if(row$status=='complete') {
        depth <- ape::node.depth.edgelength(tree)
        row$N <- ape::Ntip(tree)
        row$crown.age <- max(depth)
        cut <- row$crown.age-recent.time
        a <- depth[tree$edge[,1]]
        b <- depth[tree$edge[,2]]
        nodes <- depth[(row$N+1):length(depth)]
        row$crossing <- if(cut>0) sum(a<cut & b>=cut) else 0L
        row$old.nodes <- sum(nodes<cut)
        row$recent.nodes <- sum(nodes>=cut)
        row$old.length <- sum(pmax(0,pmin(b,cut)-a))
        row$recent.length <- sum(pmax(0,b-pmax(a,cut)))
        if(cut>0) stopifnot(row$crossing==row$old.nodes+1,
          row$recent.nodes==row$N-row$crossing)
        row$accepted <- row$N>=50 && row$crossing>=10 &&
          row$old.nodes>=10 && row$recent.nodes>=10
      }
      rows[[length(rows)+1]] <- row
    }
    attempts <- do.call(rbind, rows)
    write.csv(attempts, file.path(out,'attempts.csv'), row.names=FALSE)
    cell <- attempts[attempts$scenario==design$scenario_id[s] & attempts$boundary.fraction==fractions[f],]
    message(design$scenario_id[s], ' boundary=', fractions[f], ': ',
      sum(cell$accepted), '/', nrow(cell), ' accepted; ', sum(cell$status=='taxon_cap'), ' capped')
  }
}

# Summaries describe accepted trees; capped runs have no known final tip count. ----
summary <- list()
for(s in design$scenario_id) for(f in fractions) {
  cell <- attempts[attempts$scenario==s & attempts$boundary.fraction==f,]
  kept <- cell[cell$accepted,]
  row <- data.frame(scenario=s, boundary.fraction=f, attempted=nrow(cell),
    accepted=nrow(kept), capped=sum(cell$status=='taxon_cap'), no.tree=sum(cell$status=='no_tree'))
  for(variable in c('N','crossing','old.nodes','recent.nodes','old.length','recent.length')) {
    values <- kept[[variable]]
    for(stat in c('median','min','max','q10','q90')) {
      value <- NA_real_
      if(length(values)) value <- switch(stat, median=median(values), min=min(values),
        max=max(values), q10=unname(quantile(values,.1)), q90=unname(quantile(values,.9)))
      row[[paste(variable,stat,sep='.')]] <- value
    }
  }
  summary[[length(summary)+1]] <- row
}
summary <- do.call(rbind, summary)
write.csv(summary, file.path(out,'summary.csv'), row.names=FALSE)

library(ggplot2)
kept <- attempts[attempts$accepted,]
p <- ggplot(kept,aes(factor(boundary.fraction),N))+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=.12,height=0,seed=1),alpha=.5)+
  scale_y_log10()+facet_wrap(~scenario)+theme_classic()+
  labs(x='Boundary age / 36-unit simulation duration',y='Accepted extant tips (log scale)',
    caption='50 attempts per cell. Conditional on information rules AND completion before the 4000-taxon cap.')
ggsave(file.path(out,'tip-distributions.png'),p,width=10,height=6,dpi=160)
ggsave(file.path(out,'tip-distributions.pdf'),p,width=10,height=6)
p <- ggplot(kept,aes(factor(boundary.fraction),crossing))+
  geom_boxplot(outlier.shape=NA)+geom_point(position=position_jitter(width=.12,height=0,seed=1),alpha=.5)+
  facet_wrap(~scenario)+theme_classic()+
  labs(x='Boundary age / 36-unit simulation duration',y='Reconstructed crossing lineages')
ggsave(file.path(out,'crossing-distributions.png'),p,width=10,height=6,dpi=160)
ggsave(file.path(out,'crossing-distributions.pdf'),p,width=10,height=6)
writeLines(c('# Tree-information screening','',
  '50 fresh attempts per cell: five rate scenarios times three boundary positions. Initial state sampled uniformly. No inference.',
  'Total stem simulation duration is 36; boundary ages are 9, 18, 27 before present. These are NOT fixed fractions of reconstructed crown age.',
  'Acceptance: at least 50 extant tips, 10 reconstructed crossing lineages, and 10 reconstructed branching events per epoch.',
  'For a binary reconstructed tree spanning the boundary, old.nodes = crossing - 1 and recent.nodes = tips - crossing. Therefore the old-node rule effectively requires at least 11 crossing lineages.',
  'All distributions below are conditional on acceptance and avoiding the 4000-taxon safety cap. Capped runs are not small trees or extinction failures: they are censored and their final tip counts are unknown.',
  'Ranges and quantiles describe this screening sample, not confidence intervals. Small or zero retained samples cannot characterize a distribution.',
  '', '~~~',capture.output(print(summary[,c('scenario','boundary.fraction','accepted','capped','N.median','N.min','N.max','crossing.median','crossing.min','crossing.max')],row.names=FALSE)),'~~~'),
  file.path(out,'report.md'))
print(summary[,c('scenario','boundary.fraction','accepted','capped','N.median','N.min','N.max','crossing.median','crossing.min','crossing.max')],row.names=FALSE)
