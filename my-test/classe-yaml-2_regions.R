library(yaml)
library(diversitree)
source('R/utils-yaml.R')
source('R/utils-sse.R')

par.categories <- read_yaml_pars(yaml_file="my-test/data/2-regions.yml")
view_pars_as_array(par.categories)
formula <- make_constraints_sse(par.categories)

phy <- readRDS('my-test/data/phy-50tips.rds')
plot(phy)
lik.c <- make.classe(phy, phy$tip.state+1, 3, control=list(backend = "gslode"))
#lik.c <- make.classe(phy, phy$tip.state+1, 3)

lik.const <- constrain(lik.c, formulae = formula)
arg.const <- argnames(lik.const)
arg.const
starting.point <- starting.point.classe(phy, 3)
starting.point <-starting.point[arg.const]

mle <- find.mle(lik.const, starting.point/2, condition.surv=TRUE, keep.func=F)
mle$lnLik
# -64.70666

# #--------------------
# reg <- c("A", "B", "A.B")
# pars <- diversitree:::default.argnames.classe(3)
# args <- pars_to_arrays(pars, 3, reg)
# args
#
# par.vec <- yaml2vec("my-test/data/2-regions.yml", args)
# pars.formula <- pars2formula(par.vec$pars, pars)
# #-----

# ClaSSE td
par.categories.td <- read_yaml_pars_td("my-test/data/2-regions-td.yml")
formula.td <- make_constraints_sse_td(par.categories.td)

phy <- readRDS('my-test/data/phy-50tips.rds')
plot(phy)
#lik.c <- make.classe(phy, phy$tip.state+1, 3, control=list(backend = "gslode"))
lik.td <-make.classe.td(phy, phy$tip.state+1, k=3, n.epoch=2, control=list(backend = "gslode"))

lik.const.td <- constrain(lik.td, formulae = formula.td)
arg.const.td <- argnames(lik.const.td)
arg.const.td
#starting.point <- starting.point.classe(phy, 3)
names(starting.point) <- arg.const.td

mle.td <- find.mle(lik.const.td, starting.point/2, condition.surv=TRUE, keep.func=F)
mle.td$lnLik
mle.td$par.full
# -64.70666

lik.td(mle.td$par.full)

#--------
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

