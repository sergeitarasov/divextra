# Calibrated 5/25/50-percent pilot: 50--1000 tips. ----
library(divextra)
source('inst/benchmarks/archive/pilots/geosse-parameter-recovery/R/benchmark-functions.R')
out <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/epoch-grid-pilot'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
rates <- read.csv('inst/benchmarks/archive/pilots/geosse-parameter-recovery/config/scenarios-linked.csv')
fractions <- c(.05,.25,.50)
durations <- c(18,24,30,36,42,48)

# Reused for independent calibration and final simulation samples. ----
simulate.one <- function(row,total,fraction,seed) {
  set.seed(seed)
  root <- sample(1:3,1)
  pars <- make_simulation_parameters(row)
  recent <- total*fraction
  old <- total-recent
  tb <- make.tree.classe.td(pars$classe,k=3,max.taxa=1000,
    max.t1=old,max.t2=total,x0=root,single.lineage=TRUE)
  status <- 'complete'
  if(is.null(tb)) status <- 'no_tree'
  if(status=='complete' && !any(!tb$split & !tb$extinct)) status <- 'no_tree'
  if(status=='complete' && (abs(attr(tb,'t.epoch1')-old)>1e-8 ||
     abs(attr(tb,'t.total')-total)>1e-8)) status <- 'taxon_cap'
  if(status=='complete') {
    tree <- tryCatch(table2tree(tb),error=function(e) NULL)
    if(is.null(tree) || ape::Ntip(tree)<2) status <- 'no_tree'
  }
  diag <- data.frame(N=NA_integer_,boundary.lineages=NA_integer_,
    old.nodes=NA_integer_,recent.nodes=NA_integer_,old.length=NA_real_,
    recent.length=NA_real_,crown.age=NA_real_,boundary.age=recent)
  if(status=='complete') {
    depth <- ape::node.depth.edgelength(tree)
    cut <- max(depth)-recent
    a <- depth[tree$edge[,1]]
    b <- depth[tree$edge[,2]]
    nodes <- depth[(ape::Ntip(tree)+1):length(depth)]
    diag$N <- ape::Ntip(tree)
    diag$crown.age <- max(depth)
    diag$boundary.lineages <- if(cut>0) sum(a<cut & b>=cut) else 0L
    diag$old.nodes <- sum(nodes<cut)
    diag$recent.nodes <- sum(nodes>=cut)
    diag$old.length <- sum(pmax(0,pmin(b,cut)-a))
    diag$recent.length <- sum(pmax(0,b-pmax(a,cut)))
    if(cut>0) stopifnot(diag$old.nodes==diag$boundary.lineages-1,
      diag$recent.nodes==diag$N-diag$boundary.lineages)
    if(diag$N<50 || diag$N>1000 || diag$boundary.lineages<10 ||
       diag$old.nodes<10 || diag$recent.nodes<10) status <- 'information_filter'
  }
  result <- list(status=status,diagnostics=diag)
  if(status!='complete') return(result)
  hist <- attr(tb,'hist')
  end <- tb$start+tb$len
  for(epoch in c('old','recent')) {
    in.epoch <- if(epoch=='old') end<old else end>=old
    h.epoch <- if(epoch=='old') hist$t<old else hist$t>=old
    diag[[paste0(epoch,'.sim.speciation')]] <- sum(tb$split & in.epoch)
    diag[[paste0(epoch,'.sim.extinction')]] <- sum(tb$extinct & in.epoch)
    diag[[paste0(epoch,'.sim.dispersal')]] <- sum(h.epoch & hist$to==3)
    diag[[paste0(epoch,'.sim.contraction')]] <- sum(h.epoch & hist$from==3)
  }
  diag$n.A <- sum(tree$tip.state==1)
  diag$n.B <- sum(tree$tip.state==2)
  diag$n.AB <- sum(tree$tip.state==3)
  result$diagnostics <- diag
  result$tree <- tree
  result$table <- tb
  result$root <- root
  result
}

# Calibration uses no inference and no recovery-based selection. ----
if(!file.exists(file.path(out,'design.csv'))) {
  calibration <- list()
  design <- list()
  cell <- 0
  for(s in seq_len(nrow(rates))) for(f in seq_along(fractions)) {
    cell <- cell+1
    scores <- numeric(length(durations))
    medians <- rep(Inf,length(durations))
    for(t in seq_along(durations)) {
      sizes <- numeric()
      for(attempt in 1:30) {
        seed <- 16000000+cell*10000+t*100+attempt
        z <- simulate.one(rates[s,],durations[t],fractions[f],seed)
        calibration[[length(calibration)+1]] <- data.frame(
          scenario=rates$scenario_id[s],fraction=fractions[f],duration=durations[t],
          seed=seed,status=z$status,N=z$diagnostics$N,
          crossing=z$diagnostics$boundary.lineages)
        if(z$status=='complete') sizes <- c(sizes,z$diagnostics$N)
      }
      scores[t] <- length(sizes)
      if(length(sizes)) medians[t] <- median(sizes)
    }
    # Prefer acceptance probability; break ties by proximity to 200 tips.
    chosen <- order(-scores,abs(log(medians/200)))[1]
    row <- rates[s,]
    row$rate.scenario <- row$scenario_id
    row$scenario_id <- paste0(row$scenario_id,'_b',sprintf('%02d',100*fractions[f]))
    row$fraction <- fractions[f]
    row$total.time <- durations[chosen]
    row$old.time <- durations[chosen]*(1-fractions[f])
    row$recent.time <- durations[chosen]*fractions[f]
    row$calibration.accepted <- scores[chosen]
    design[[cell]] <- row
    write.csv(do.call(rbind,calibration),file.path(out,'calibration-attempts.csv'),row.names=FALSE)
    message('Calibrated ',row$scenario_id,': duration=',row$total.time,
      '; accepted=',scores[chosen],'/30')
  }
  design <- do.call(rbind,design)
  write.csv(design,file.path(out,'design.csv'),row.names=FALSE)
} else design <- read.csv(file.path(out,'design.csv'))

# Five fresh qualifying trees per cell; all rejections are retained. ----
if(!file.exists(file.path(out,'datasets.rds'))) {
  datasets <- list()
  attempts <- list()
  for(cell in seq_len(nrow(design))) {
    row <- design[cell,]
    retained <- 0
    for(attempt in 1:2000) {
      seed <- 18000000+cell*10000+attempt
      z <- simulate.one(row,row$total.time,row$fraction,seed)
      attempts[[length(attempts)+1]] <- data.frame(scenario=row$scenario_id,
        seed=seed,status=z$status,N=z$diagnostics$N,crossing=z$diagnostics$boundary.lineages)
      if(z$status!='complete') next
      retained <- retained+1
      datasets[[length(datasets)+1]] <- list(id=paste0(row$scenario_id,'-',retained),
        scenario=row$scenario_id,rate.scenario=row$rate.scenario,fraction=row$fraction,
        links=c(s=row$link_s,d=row$link_d,x=row$link_x),replicate=retained,
        tree=z$tree,table=z$table,seed=seed,root=z$root,
        truth=c(s_recent=row$s_recent,d_recent=row$d_recent,x_recent=row$x_recent,
          s_old=row$s_old,d_old=row$d_old,x_old=row$x_old),
        recent.time=row$recent.time,diagnostics=z$diagnostics)
      if(retained==5) break
    }
    write.csv(do.call(rbind,attempts),file.path(out,'simulation-attempts.csv'),row.names=FALSE)
    saveRDS(datasets,file.path(out,'datasets-partial.rds'))
    message('Generated ',row$scenario_id,': ',retained,'/5 in ',attempt,' attempts')
    if(retained<5) stop('Cell incomplete; calibration needs review: ',row$scenario_id)
  }
  saveRDS(datasets,file.path(out,'datasets.rds'))
} else datasets <- readRDS(file.path(out,'datasets.rds'))
stopifnot(length(datasets)==75)

# Common parameter order: s, d, x (recent then old for episodic fits). ----
fit.one <- function(z,model,start) {
  library(divextra)
  if(model=='geosse') {
    states <- setNames(c(1L,2L,0L)[z$tree$tip.state],names(z$tree$tip.state))
    raw <- diversitree::make.geosse(z$tree,states,strict=FALSE,
      sampling.f=c(1,1,1),control=list(tol=1e-10))
    evaluate <- function(p) raw(c(p[1],p[1],p[1],p[3],p[3],p[2],p[2]),
      root=diversitree::ROOT.OBS,condition.surv=TRUE)
  } else {
    file <- 'model-episodic.yml'
    cfg <- read_yaml_pars_td(file.path('inst/benchmarks/archive/pilots/geosse-parameter-recovery/config',file))
    cfg$epoch.times <- z$recent.time
    raw <- make.classe.td(z$tree,z$tree$tip.state,k=3,n.epoch=2,strict=FALSE,
      sampling.f=c(1,1,1),control=list(tol=1e-10))
    raw <- diversitree::constrain(raw,formulae=make_constraints_sse_td(cfg))
    expected <- c('lambda111.1','mu1.1','q13.1','lambda111.2','mu1.2','q13.2')
    stopifnot(identical(diversitree::argnames(raw),expected))
    evaluate <- function(p) {
      full <- setNames(numeric(6),names(z$truth))
      for(rate in c('s','d','x')) {
        full[paste0(rate,'_recent')] <- p[paste0(rate,'_recent')]
        old.name <- if(z$links[rate]) paste0(rate,'_recent') else paste0(rate,'_old')
        full[paste0(rate,'_old')] <- p[old.name]
      }
      values <- full[c(1,3,2,4,6,5)]
      raw(values,root=diversitree::ROOT.OBS,condition.surv=TRUE)
    }
  }
  parameter.names <- c('s_recent','d_recent','x_recent')
  for(rate in c('s','d','x')) if(!z$links[rate]) parameter.names <- c(parameter.names,paste0(rate,'_old'))
  if(model=='geosse') parameter.names <- c('s','d','x')
  # Distinct starting rates, no use of generating rates for optimization.
  set.seed(z$seed+start*100000)
  initial <- rep(c(.08,.08,.03),length.out=length(parameter.names))*c(.5,1,2)[start]
  initial <- initial*exp(runif(length(initial),-log(2),log(2)))
  names(initial) <- parameter.names
  objective <- function(logp) {
    value <- tryCatch(suppressWarnings(evaluate(exp(logp))),error=function(e) -Inf)
    if(!is.finite(value)) return(1e50)
    -value
  }
  begin <- proc.time()['elapsed']
  fit <- nlminb(log(initial),objective,lower=log(1e-10),upper=log(100),
    control=list(iter.max=1000,eval.max=3000,rel.tol=1e-9))
  pars <- setNames(exp(fit$par),parameter.names)
  truth <- if(model=='geosse') unname(z$truth[1:3]) else z$truth[parameter.names]
  list(id=z$id,scenario=z$scenario,model=model,start=start,initial=initial,
    par=pars,lnLik=-fit$objective,convergence=fit$convergence,message=fit$message,
    bound.hit=any(pars>99),elapsed=unname(proc.time()['elapsed']-begin),
    truth.lnLik=evaluate(truth),optimizer=fit)
}
jobs <- list()
for(i in seq_along(datasets)) {
  models <- 'linked_classe'
  if(datasets[[i]]$rate.scenario=='constant') models <- c(models,'geosse')
  for(model in models) for(start in 1:3)
    jobs[[length(jobs)+1]] <- list(dataset=i,model=model,start=start)
}
fits <- vector('list',length(jobs))
fit.file <- file.path(out,'fits.rds')
if(file.exists(fit.file)) fits <- readRDS(fit.file)
pending <- which(vapply(fits,is.null,logical(1)))
if(length(pending)) for(first in seq(1,length(pending),by=3)) {
  indices <- pending[first:min(first+2,length(pending))]
  processes <- list()
  for(j in indices) {
    job <- jobs[[j]]
    processes[[as.character(j)]] <- callr::r_bg(fit.one,
      list(z=datasets[[job$dataset]],model=job$model,start=job$start))
  }
  begin <- Sys.time()
  repeat {
    alive <- FALSE
    for(j in indices) {
      p <- processes[[as.character(j)]]
      if(p$is_alive()) {
        if(as.numeric(difftime(Sys.time(),begin,units='secs'))>180) p$kill() else alive <- TRUE
      }
    }
    if(!alive) break
    Sys.sleep(.2)
  }
  for(j in indices) {
    job <- jobs[[j]]
    z <- datasets[[job$dataset]]
    fits[[j]] <- tryCatch(processes[[as.character(j)]]$get_result(),error=function(e)
      list(id=z$id,scenario=z$scenario,model=job$model,start=job$start,error=conditionMessage(e)))
    message(j,'/',length(jobs),' ',z$id,' ',job$model,' start ',job$start,
      if(is.null(fits[[j]]$error)) paste0(' code=',fits[[j]]$convergence) else ' failed')
  }
  saveRDS(fits,fit.file)
}
saveRDS(jobs,file.path(out,'jobs.rds'))


