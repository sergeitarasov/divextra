# Independent mathematical audit of extinct side trees ----
# Run from anywhere after loading the current divextra package.
# No fitted parameters or empirical data are needed. No production code is changed.
if (!"divextra" %in% loadedNamespaces()) library(divextra)
cache.fun <- divextra:::.classe_td_extinction_cache
schedule.fun <- divextra:::.classe_td_schedule
draw.fun <- divextra:::.classe_td_draw_extinct_event

# Record a complete side tree, independent of the package graph builder ----
draw.side <- function(schedule, cache, age, state=1L, oracle=NULL) {
  stack <- list(c(age,state))
  segments <- list()
  deaths <- numeric()
  death.states <- integer()
  births <- transitions <- 0L
  while (length(stack)) {
    lineage <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    t <- lineage[1]; i <- as.integer(lineage[2])
    repeat {
      if (is.null(oracle)) {
        event <- draw.fun(t,i,numeric(),schedule,cache)
        u <- event$age
        type <- event$type
        pair <- event$pair
        to <- event$to
      } else {
        # Independent ordinary Gillespie process with epoch crossing.
        epoch <- 1L + sum(t > oracle$boundaries)
        younger <- c(0,oracle$boundaries[oracle$boundaries < t])
        bound <- max(younger)
        rr <- oracle$rates[[epoch]]
        b <- rr$birth[[i]]
        dest <- setdiff(seq_along(rr$mu),i)
        w <- c(b[,3],rr$mu[i],rr$Q[i,dest])
        dt <- rexp(1,sum(w))
        if (dt >= t-bound) {
          segments[[length(segments)+1L]] <- c(start=t,end=bound,state=i)
          if (bound == 0) return(NULL) # Reject at the first surviving descendant.
          t <- bound
          next
        }
        u <- t-dt
        z <- sample.int(length(w),1,prob=w)
        type <- if(z <= nrow(b)) "speciation" else if(z==nrow(b)+1L) "extinction" else "transition"
        pair <- if(type=="speciation") b[z,1:2] else NULL
        to <- if(type=="transition") dest[z-nrow(b)-1L] else NULL
      }
      segments[[length(segments)+1L]] <- c(start=t,end=u,state=i)
      if (type == "extinction") {
        deaths <- c(deaths,u); death.states <- c(death.states,i)
        break
      }
      if (type == "speciation") {
        births <- births+1L
        stack[[length(stack)+1L]] <- c(u,pair[1])
        stack[[length(stack)+1L]] <- c(u,pair[2])
        break
      }
      transitions <- transitions+1L
      t <- u; i <- to
    }
  }
  segments <- do.call(rbind,segments)
  stopifnot(length(deaths)==births+1L,all(deaths>0),
            all(segments[,"start"]>segments[,"end"]))
  alive <- segments[,"start"] > age/2 & segments[,"end"] <= age/2
  c(births=births, tips=length(deaths), transitions=transitions,
    length=sum(segments[,"start"]-segments[,"end"]),
    midpoint=sum(alive), midpoint.state2=sum(alive & segments[,"state"]==2),
    death.state2=sum(death.states==2),
    extinction.elapsed=age-min(deaths), no.birth=as.numeric(births==0))
}

# One-state finite-time birth-death formulas ----
bd.E <- function(t,lambda,mu) {
  if(lambda==mu) return(mu*t/(1+mu*t))
  mu*(-expm1(-(lambda-mu)*t))/(lambda-mu*exp(-(lambda-mu)*t))
}
bd.mean <- function(s,T,lambda,mu) {
  p0 <- bd.E(s,lambda,mu)
  beta <- lambda/mu*p0
  z <- bd.E(T-s,lambda,mu)
  (1-p0)*(1-beta)*z/(1-beta*z)^2/bd.E(T,lambda,mu)
}
summarise.oracle <- function(draws, expected, model) {
  do.call(rbind,lapply(names(expected),function(metric) {
    v <- draws[metric,]
    se <- sd(v)/sqrt(length(v))
    data.frame(model=model,metric=metric,estimated=mean(v),expected=expected[[metric]],
               MC.SE=se,z=if(se>0) (mean(v)-expected[[metric]])/se else NA_real_)
  }))
}
set.seed(260923)
N <- 2000L
analytic.results <- list()
for (case in list(c(lambda=0,mu=.7),c(lambda=.3,mu=.7),
                  c(lambda=.7,mu=.3),c(lambda=.5,mu=.5))) {
  lambda <- case[[1]]; mu <- case[[2]]; T <- 3
  schedule <- schedule.fun(c(lambda111.1=lambda,mu1.1=mu),1,"A")
  cache <- cache.fun(schedule,T)
  draws <- replicate(N,draw.side(schedule,cache,T))
  mean.length <- integrate(function(s) bd.mean(s,T,lambda,mu),0,T)$value
  mean.births <- integrate(function(s)
    lambda*bd.E(T-s,lambda,mu)*bd.mean(s,T,lambda,mu),0,T)$value
  expected <- c(tips=1+mean.births,length=mean.length,
    midpoint=bd.mean(T/2,T,lambda,mu),
    no.birth=mu/(lambda+mu)*(-expm1(-(lambda+mu)*T))/bd.E(T,lambda,mu))
  name <- paste0("lambda=",lambda,",mu=",mu)
  analytic.results[[name]] <- summarise.oracle(draws,expected,name)
  # CDF of total clade extinction time, conditional on extinction by T.
  pit <- bd.E(draws["extinction.elapsed",],lambda,mu)/bd.E(T,lambda,mu)
  distance <- max(abs(sort(pit)-(seq_len(N)-.5)/N))
  cat(name,"extinction-time uniform ECDF distance",distance,"\n")
}
analytic.results <- do.call(rbind,analytic.results)
print(analytic.results,row.names=FALSE)

# Independent multistate, two-epoch forward-rejection oracle ----
young <- list(mu=c(.7,.5),Q=matrix(c(0,.25,.1,0),2,byrow=TRUE),
  birth=list(cbind(c(1,1,2),c(1,2,2),c(.2,.15,.05)),
             cbind(c(1,1,2),c(1,2,2),c(.03,.2,.1))))
old <- young
old$mu <- c(1.05,.35)
old$Q <- matrix(c(0,.1,.3,0),2,byrow=TRUE)
old$birth <- lapply(old$birth,function(b){b[,3]<-b[,3]*.8;b})
oracle <- list(boundaries=1,rates=list(young,old))
pack <- function(r,e) {
  p <- c(lambda111=r$birth[[1]][1,3],lambda112=r$birth[[1]][2,3],
    lambda122=r$birth[[1]][3,3],lambda211=r$birth[[2]][1,3],
    lambda212=r$birth[[2]][2,3],lambda222=r$birth[[2]][3,3],
    mu1=r$mu[1],mu2=r$mu[2],q12=r$Q[1,2],q21=r$Q[2,1])
  setNames(p,paste0(names(p),".",e))
}
schedule <- schedule.fun(c(pack(young,1),pack(old,2),t.1=1),2,c("A","B"))
cache <- cache.fun(schedule,3)
direct <- replicate(N,draw.side(schedule,cache,3))
reference <- matrix(NA_real_,nrow(direct),N,dimnames=list(rownames(direct),NULL))
attempts <- 0L
for(i in seq_len(N)) repeat {
  attempts <- attempts+1L
  x <- draw.side(NULL,NULL,3,oracle=oracle)
  if(!is.null(x)) {reference[,i]<-x;break}
}
multistate.results <- do.call(rbind,lapply(rownames(direct),function(metric) {
  a<-direct[metric,];b<-reference[metric,]
  se<-sqrt(var(a)/N+var(b)/N)
  data.frame(metric=metric,direct=mean(a),rejection=mean(b),MC.SE=se,
             z=if(se>0)(mean(a)-mean(b))/se else NA_real_)
}))
cat("Multistate oracle attempts",attempts,"for",N,"accepted trees\n")
print(multistate.results,row.names=FALSE)

# Public backbone: expected hidden-birth count under OUR augmentation ----
for (mu in c(1,.05)) {
maps <- make.simmap.classe.td(ape::read.tree(text="(a:2,b:2);"),
  c(lambda111.1=.5,mu1.1=mu),tip.weights=c(a=1,b=1),
  node.probs=matrix(1,1,1,dimnames=list("A","3")),k=1,state.labels="A",
  side.lineage.mode="extinct_clades",nsim=1000,seed=120)
births <- vapply(maps,function(s)sum(s$events$event=="hidden_speciation"),0)
cat("Backbone hidden births: mean",mean(births),"expected",.5*4,
    "MC SE",sd(births)/sqrt(length(births)),"\n")
cat("mu",mu,"joint tree-conditioned expected hidden births (different target):",
    2*integrate(function(t)2*.5*bd.E(t,.5,mu),0,2)$value,"\n")
}

# Stress numerical interpolation separately from Monte Carlo ----
stress.results <- do.call(rbind,lapply(c(.7,100,10000),function(mu) {
  schedule<-schedule.fun(c(lambda111.1=0,mu1.1=mu),1,"A")
  cache<-tryCatch(cache.fun(schedule,3),error=function(e)e)
  if(inherits(cache,"error")) return(data.frame(mu=mu,max.E.error=NA_real_,
    status=conditionMessage(cache)))
  times<-sort(unique(c(seq(0,3,length.out=10001),seq(0,5/mu,length.out=10001))))
  times<-times[times<=3]
  data.frame(mu=mu,max.E.error=max(abs(vapply(times,function(t)cache$E(t)[1],0)-
    (-expm1(-mu*times)))),status="completed")
}))
print(stress.results,row.names=FALSE)

# A pre-recorded follow-up to the initial critical-case discrepancy ----
# Retain the initial 2,000-draw result above; use an independent seed here.
schedule <- schedule.fun(c(lambda111.1=.5,mu1.1=.5),1,"A")
cache <- cache.fun(schedule,3)
set.seed(73480)
critical.draws <- replicate(20000,draw.side(schedule,cache,3))
critical.results <- summarise.oracle(critical.draws,
  c(tips=1.3,length=1.2,midpoint=.35,no.birth=.5*(-expm1(-3))/.6),
  "critical independent 20000")
print(critical.results,row.names=FALSE)
