# Prototype benchmark: exact one-state, extinction-only backbone ----
# Standalone base R. This does NOT modify or call the production backbone sampler.
# Each benchmark is one stem of duration 2 with exactly one surviving descendant.
# Two such stems form the two-tip example discussed in the audit.
extinction <- function(t,lambda,mu) {
  if(lambda==mu) return(lambda*t/(1+lambda*t))
  mu*(-expm1(-(lambda-mu)*t))/(lambda-mu*exp(-(lambda-mu)*t))
}
cumulative <- function(t,lambda,mu) {
  if(lambda==mu) return(2*(lambda*t-log1p(lambda*t)))
  r<-lambda-mu
  2*lambda*t-2*log((lambda*exp(r*t)-mu)/r)
}

# Proposed conditional-backbone rule: a nonhomogeneous Poisson process ----
draw.conditioned <- function(T,lambda,mu) {
  total <- cumulative(T,lambda,mu)
  n <- rpois(1,total)
  if(!n) return(numeric())
  targets <- runif(n,0,total)
  sort(vapply(targets,function(z)
    uniroot(function(t)cumulative(t,lambda,mu)-z,c(0,T),tol=1e-10)$root,0))
}

# Independent unconditioned forward birth-death process ----
# Record the hidden births on the path to the sole survivor, if there is one.
draw.forward <- function(age,lambda,mu) {
  dt <- rexp(1,lambda+mu)
  if(dt >= age) return(list(n=1L,ages=numeric()))
  if(runif(1) < mu/(lambda+mu)) return(list(n=0L,ages=numeric()))
  age <- age-dt
  left <- draw.forward(age,lambda,mu)
  if(left$n >= 2L) return(list(n=2L,ages=numeric()))
  right <- draw.forward(age,lambda,mu)
  n <- left$n+right$n
  if(n==1L) return(list(n=1L,ages=c(age,if(left$n==1L)left$ages else right$ages)))
  list(n=min(n,2L),ages=numeric())
}

# Compare full count distributions and conditional event ages ----
benchmarks <- lapply(c(.05,1),function(mu) {
  lambda <- .5; T <- 2; N <- 5000L
  set.seed(1401)
  direct.time <- system.time(direct<-replicate(N,draw.conditioned(T,lambda,mu),simplify=FALSE))["elapsed"]
  set.seed(2802)
  attempts <- 0L
  reference <- vector("list",N)
  reference.time <- system.time(for(i in seq_len(N)) repeat {
    attempts <- attempts+1L
    x <- draw.forward(T,lambda,mu)
    if(x$n==1L){reference[[i]]<-x$ages;break}
  })["elapsed"]
  a <- lengths(direct); b <- lengths(reference)
  expected <- cumulative(T,lambda,mu)
  se <- sqrt(var(a)/N+var(b)/N)
  ages <- unlist(direct)
  ref.ages <- unlist(reference)
  ecdf.distance <- function(ages) {
    if(!length(ages))return(NA_real_)
    u<-sort(cumulative(ages,lambda,mu)/expected)
    max(abs(u-(seq_along(u)-.5)/length(u)))
  }
  count.distance <- function(x)max(abs(vapply(0:8,function(n)mean(x<=n),0)-ppois(0:8,expected)))
  # An unconditioned process ends in exactly one survivor with this probability.
  e <- extinction(T,lambda,mu)
  p1 <- (1-e)*(1-lambda/mu*e)
  stopifnot(abs(mean(a)-expected)<5*sqrt(expected/N),
            abs(mean(b)-expected)<5*sqrt(expected/N),
            abs(mean(a)-mean(b))<5*se)
  data.frame(lambda=lambda,mu=mu,n=N,expected=expected,
    direct.mean=mean(a),forward.mean=mean(b),comparison.z=(mean(a)-mean(b))/se,
    direct.count.CDF.error=count.distance(a),forward.count.CDF.error=count.distance(b),
    direct.age.CDF.error=ecdf.distance(ages),forward.age.CDF.error=ecdf.distance(ref.ages),
    forward.attempts=attempts,expected.acceptance=p1,observed.acceptance=N/attempts,
    direct.seconds=unname(direct.time),forward.seconds=unname(reference.time))
})
results <- do.call(rbind,benchmarks)
print(results,row.names=FALSE)

# Two epochs: independent forward process versus conditional thinning ----
lambda.epoch <- c(.2,.7)
mu.epoch <- c(.8,.1)
boundary <- 1
T <- 3
epoch.E <- function(age) {
  young <- extinction(pmin(age,boundary),lambda.epoch[1],mu.epoch[1])
  z <- extinction(boundary,lambda.epoch[1],mu.epoch[1])
  p0 <- extinction(pmax(age-boundary,0),lambda.epoch[2],mu.epoch[2])
  beta <- lambda.epoch[2]/mu.epoch[2]*p0
  # Compose the older birth-death generating function with younger extinction.
  old <- p0+(1-p0)*(1-beta)*z/(1-beta*z)
  ifelse(age<=boundary,young,old)
}
epoch.rate <- function(age)2*lambda.epoch[1+(age>boundary)]*epoch.E(age)
epoch.cumulative <- function(age)vapply(age,function(t) {
  integrate(epoch.rate,0,min(t,boundary),rel.tol=1e-11)$value+
    if(t>boundary)integrate(epoch.rate,boundary,t,rel.tol=1e-11)$value else 0
},0)
draw.epoch.conditional <- function() {
  # Exact Poisson thinning: E<=1 gives a constant dominating intensity.
  ceiling.rate <- 2*max(lambda.epoch)
  ages <- runif(rpois(1,ceiling.rate*T),0,T)
  ages[runif(length(ages))<epoch.rate(ages)/ceiling.rate]
}
draw.epoch.forward <- function(age) {
  repeat {
    epoch <- 1+(age>boundary)
    limit <- if(epoch==2)boundary else 0
    dt <- rexp(1,lambda.epoch[epoch]+mu.epoch[epoch])
    if(dt>=age-limit) {
      if(limit==0)return(list(n=1L,ages=numeric()))
      age<-limit
      next
    }
    age<-age-dt
    if(runif(1)<mu.epoch[epoch]/(lambda.epoch[epoch]+mu.epoch[epoch]))
      return(list(n=0L,ages=numeric()))
    left<-draw.epoch.forward(age)
    if(left$n>=2L)return(list(n=2L,ages=numeric()))
    right<-draw.epoch.forward(age)
    n<-left$n+right$n
    if(n==1L)return(list(n=1L,ages=c(age,if(left$n==1L)left$ages else right$ages)))
    return(list(n=min(n,2L),ages=numeric()))
  }
}
set.seed(7319)
N<-5000L
direct.time<-system.time(direct<-replicate(N,draw.epoch.conditional(),simplify=FALSE))["elapsed"]
set.seed(9107)
reference<-vector("list",N)
attempts<-0L
forward.time<-system.time(for(i in seq_len(N))repeat {
  attempts<-attempts+1L
  x<-draw.epoch.forward(T)
  if(x$n==1L){reference[[i]]<-x$ages;break}
})["elapsed"]
expected<-epoch.cumulative(T)
a<-lengths(direct);b<-lengths(reference)
se<-sqrt(var(a)/N+var(b)/N)
stopifnot(abs(mean(a)-expected)<5*sqrt(expected/N),
          abs(mean(b)-expected)<5*sqrt(expected/N),
          abs(mean(a)-mean(b))<5*se)
epoch.results<-data.frame(n=N,expected=expected,direct.mean=mean(a),
  forward.mean=mean(b),comparison.z=(mean(a)-mean(b))/se,
  direct.seconds=unname(direct.time),forward.seconds=unname(forward.time),
  forward.attempts=attempts)
print(epoch.results,row.names=FALSE)
ages<-unlist(direct)
ref.ages<-unlist(reference)
cat("Expected fraction of hidden births in young epoch:",epoch.cumulative(boundary)/expected,
    "direct:",mean(ages<=boundary),"forward:",mean(ref.ages<=boundary),"\n")
