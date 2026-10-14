# Fixed-time rate variants and classical GeoSSE control.
# Run from the DivExtra repository root. Existing completed fits are reused.
library(divextra)
source('inst/benchmarks/archive/pilots/geosse-parameter-recovery/R/benchmark-functions.R')
out <- 'inst/benchmarks/archive/pilots/geosse-parameter-recovery/results/linked-rate-pilot'
dir.create(out,recursive=TRUE,showWarnings=FALSE)

# Correctly specified 6/4/4/4/3-parameter models. ----
design <- read.csv('inst/benchmarks/archive/pilots/geosse-parameter-recovery/config/scenarios-linked.csv')
design$old.time <- 30
design$recent.time <- 6
stopifnot(design$n_parameters == 6-rowSums(design[c('link_s','link_d','link_x')]))
write.csv(design,file.path(out,'design.csv'),row.names=FALSE)

datasets <- list()
attempts <- list()
for(s in seq_len(nrow(design))) {
  row <- design[s,]
  parameters <- make_simulation_parameters(row)
  retained <- 0
  for(attempt in 1:500) {
    seed <- 12700000+s*10000+attempt
    set.seed(seed)
    root <- sample(1:3,1)
    tb <- make.tree.classe.td(parameters$classe,k=3,max.taxa=4000,
      max.t1=row$old.time,max.t2=row$old.time+row$recent.time,x0=root,single.lineage=TRUE)
    status <- 'accepted'
    tree <- NULL
    if(is.null(tb)) status <- 'no_tree'
    if(status=='accepted' && (abs(attr(tb,'t.epoch1')-row$old.time)>1e-8 ||
      abs(attr(tb,'t.total')-row$old.time-row$recent.time)>1e-8)) status <- 'taxon_limit'
    if(status=='accepted') {
      tree <- tryCatch(table2tree(tb),error=function(e) NULL)
      if(is.null(tree) || ape::Ntip(tree)<2) status <- 'no_tree'
    }
    if(status=='accepted') {
      depth <- ape::node.depth.edgelength(tree)
      cut <- max(depth)-row$recent.time
      if(cut<=0) status <- 'no_old_crown'
    }
    attempts[[length(attempts)+1]] <- data.frame(scenario=row$scenario_id,attempt=attempt,status=status)
    if(status!='accepted') next
    retained <- retained+1
    a <- depth[tree$edge[,1]]
    b <- depth[tree$edge[,2]]
    datasets[[length(datasets)+1]] <- list(id=paste0(row$scenario_id,'-',retained),
      scenario=row$scenario_id,links=c(s=row$link_s,d=row$link_d,x=row$link_x),replicate=retained,tree=tree,table=tb,seed=seed,root=root,
      truth=c(s_recent=row$s_recent,d_recent=row$d_recent,x_recent=row$x_recent,
        s_old=row$s_old,d_old=row$d_old,x_old=row$x_old),recent.time=row$recent.time,
      diagnostics=data.frame(N=ape::Ntip(tree),boundary.lineages=sum(a<cut & b>=cut),
        old.length=sum(pmax(0,pmin(b,cut)-a)),
        recent.length=sum(pmax(0,b-pmax(a,cut))),
        crown.age=max(depth),boundary.age=row$recent.time,
        n.A=sum(tree$tip.state==1),n.B=sum(tree$tip.state==2),n.AB=sum(tree$tip.state==3),
        old.nodes=sum(depth[(ape::Ntip(tree)+1):length(depth)]<cut),
        recent.nodes=sum(depth[(ape::Ntip(tree)+1):length(depth)]>=cut)))
    # Complete-history events; these are not all visible in the reconstructed tree.
    end <- tb$start+tb$len
    hist <- attr(tb,'hist')
    for(epoch in c('old','recent')) {
      in.epoch <- if(epoch=='old') end<row$old.time else end>=row$old.time
      h.epoch <- if(epoch=='old') hist$t<row$old.time else hist$t>=row$old.time
      z <- datasets[[length(datasets)]]
      z$diagnostics[[paste0(epoch,'.sim.speciation')]] <- sum(tb$split & in.epoch)
      z$diagnostics[[paste0(epoch,'.sim.extinction')]] <- sum(tb$extinct & in.epoch)
      z$diagnostics[[paste0(epoch,'.sim.dispersal')]] <- sum(h.epoch & hist$to==3)
      z$diagnostics[[paste0(epoch,'.sim.contraction')]] <- sum(h.epoch & hist$from==3)
      datasets[[length(datasets)]] <- z
    }
    if(retained==5) break
  }
  stopifnot(retained==5)
}
saveRDS(datasets,file.path(out,'datasets.rds'))
write.csv(do.call(rbind,attempts),file.path(out,'simulation-attempts.csv'),row.names=FALSE)

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
  if(datasets[[i]]$scenario=='constant') models <- c(models,'geosse')
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

