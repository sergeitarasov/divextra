library(yaml)
library(diversitree)
library(divextra)
source('R/utils-yaml.R')
source('R/utils-sse.R')
source('R/sim-classe-td.R')

## "sA"  "sB" "sAB" "xA"  "xB"  "dA"  "dB"
pars.ge <- matrix(
  c(0.1, 0.1, 0.1,  0, 0,  0.1, 0.1,
    0.3, 0.3, 0.3,  0, 0,  0.1, 0.1
  ),
  2, 7, byrow = TRUE)
colnames(pars.ge) <- diversitree:::default.argnames.geosse()
pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))
print(pars.cl)

set.seed(123)
tb <- make.tree.classe.td(pars.cl, k=3, max.t1=10, max.t2=15, x0=1, single.lineage=TRUE)
print(tb)
phy <- table2tree(tb)
plot(phy)

# the tree from epoch 1 only
tb.e1 <-  attr(tb, "info.epoch1")
phy.e1 <- diversitree:::me.to.ape.bisse(tb.e1[-1,], tb.e1$state[1])
phy.e1 <- diversitree::prune(phy.e1) # prune extinct
plot(phy.e1)

# the final tree (after epoch 2)
phy$t.epoch1 # time of epoch 1 (=max.t1)
phy$t.total # total time (=max.t2)
phy$epoch1.extant.tree.depth # tree depth from epoch 1
phy$epoch12.tree.depth # total tree depth
phy$epoch1.ntip # N of tips in epoch 1 only
Ntip(phy) # total N of tips

# The model changes regimes at (if counting from the tips for the infference):
phy$t.regime.change # := phy$epoch12.tree.depth - phy$epoch1.extant.tree.depth

#---- simulate GeoSSE td

  # "sA"  "sB" "sAB" "xA"  "xB"  "dA"  "dB"
pars.ge <- matrix(
  c(0.1, 0.1, 0.1,  0.001, 0.001,  0.1, 0.1,
    0.3, 0.3, 0.3,  0.001, 0.001,  0.1, 0.1
  ),
  2, 7, byrow = TRUE)
colnames(pars.ge) <- diversitree:::default.argnames.geosse()
pars.ge
pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))
pars.cl

tb <- make.tree.classe.td(pars.cl, k=3, max.t1=30, max.t2=35, x0=1, single.lineage=TRUE)
tb
phy <- table2tree(tb)
phy$epoch1.ntip
Ntip(phy)- phy$epoch1.ntip
plot(phy)
# tb.e1 <-  attr(tb, "info.epoch1")
# phy.e1 <- diversitree:::me.to.ape.bisse(tb.e1[-1,], tb.e1$state[1])
# phy.e1 <- diversitree::prune(phy.e1)

phy$t.epoch1
phy$t.total
phy$epoch1.extant.tree.depth
phy$epoch12.tree.depth
phy$epoch1.ntip
Ntip(phy)
#phy$tip.state
#length(phy$tip.state)

# !!!since in simulations epoch start from the root while in infference from tip the model swith is at:
phy$epoch12.tree.depth - phy$epoch1.extant.tree.depth


#------ Inference Classe TD
par.categories.td <- read_yaml_pars_td("my-test/data/2-regions-td-test.yml")
par.categories.td
pars.cl
formula.td <- make_constraints_sse_td(par.categories.td)

#phy <- readRDS('my-test/data/phy-50tips.rds')
plot(phy)
lik.td <-make.classe.td(phy, phy$tip.state, k=3, n.epoch=2, control=list(backend = "gslode"))

lik.const.td <- constrain(lik.td, formulae = formula.td)
arg.const.td <- argnames(lik.const.td)
arg.const.td
init <- starting.point.classe_td(phy, k=3, n.epoch=2)
starting.point <- init[arg.const.td]
starting.point

mle.td <- find.mle(lik.const.td, starting.point, condition.surv=TRUE, keep.func=F)
mle.td$lnLik
mle.td$par
get_aic(mle.td$lnLik, 4) # 3994.674
mle.td$par.full %>% round(., 3)
# lik.td(mle.td$par.full)



#------ Inference Classe  TD Constant
par.categories.td <- read_yaml_pars_td("my-test/data/2-regions-td-test-constant.yml")
par.categories.td
pars.cl
formula.td <- make_constraints_sse_td(par.categories.td)

plot(phy)
lik.td <-make.classe.td(phy, phy$tip.state, k=3, n.epoch=2, control=list(backend = "gslode"))

lik.const.td <- constrain(lik.td, formulae = formula.td)
arg.const.td <- argnames(lik.const.td)
arg.const.td
init <- starting.point.classe_td(phy, k=3, n.epoch=2)
starting.point <- init[arg.const.td]
starting.point

mle.td <- find.mle(lik.const.td, starting.point, condition.surv=TRUE, keep.func=F)
mle.td$lnLik
mle.td$par
get_aic(mle.td$lnLik, 3) # 4049.455
mle.td$par.full %>% round(., 3)
# lik.td(mle.td$par.full)

#------ Inference Geosse

# my order is (A, B, AB), geosse order is (AB, A, B)
new.state <- plyr::mapvalues(phy$tip.state, from = c(1,2,3), to=c(2,3,1))

lik.ge <- make.geosse(phy, new.state-1)
lik.d <- constrain(lik.ge, dB ~ dA, sB~sA, sAB~sA, xB~xA)
aa <- argnames(lik.d)
p <- starting.point.geosse(phy)[aa]

mle.ge <- find.mle(lik.d, p, method="subplex")
mle.ge$lnLik
mle.ge$par
get_aic(mle.ge$lnLik, 3) # 4049.455
