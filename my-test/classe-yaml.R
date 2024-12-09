library(yaml)
library(diversitree)
source('R/utils-yaml.R')
source('R/utils-sse.R')



#---
reg <- c("A", "M", "Mr", "OP", "S", "A.M", "A.Mr", "A.OP", "A.S", "Mr.S")
pars <- diversitree:::default.argnames.classe(10)
args <- pars_to_arrays(pars,10, reg)
args

par.vec <- yaml2vec("my-test/data/Afr-to-Maur.yaml", args)
pars.formula <- pars2formula(par.vec, pars)


lik.c <- make.classe(phy, dt, 10, strict=F)
lik.const <- constrain(lik.c, formulae = pars.formula)
arg.const <- argnames(lik.const)
arg.const
starting.point <- starting.point.classe(phy, 10)
starting.point <-starting.point[arg.const]

mle10 <- find.mle(lik.const, starting.point, method="subplex", keep.func=F, root=ROOT.OBS, condition.surv=TRUE)
mle10$par
mle10$lnLik
get_aic(mle10$lnLik, length(arg.const))
# 737.7395

