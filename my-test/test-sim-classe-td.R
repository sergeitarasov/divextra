library(yaml)
library(diversitree)
library(divextra)
source('R/utils-yaml.R')
source('R/utils-sse.R')
source('R/sim-classe-td.R')

#---- simulate Classe td

pars.tb <- matrix(
  c(1:27/10,
    1:27/100
  ),
  2, 27, byrow = T)
colnames(pars.tb) <- diversitree:::default.argnames.classe(3)

pars.tb


tb <- make.tree.classe.td(pars.tb, k=3, max.t1=1, max.t2=3, x0=1, single.lineage=TRUE)
tb
phy <- table2tree(tb)
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
phy$tip.state
length(phy$tip.state)


#---- simulate GeoSSE td

  # "sA"  "sB"  "sAB" "xA"  "xB"  "dA"  "dB"
pars.ge <- matrix(
  c(0.1, 0.1, 0.1,  0.001, 0.001,  0.1, 0.1,
    0.3, 0.3, 0.3,  0.001, 0.001,  0.1, 0.1
  ),
  2, 7, byrow = T)
colnames(pars.ge) <- diversitree:::default.argnames.geosse()
pars.ge
pars.cl <- t(apply(pars.ge, 1, diversitree:::pars.ge.to.cl))

tb <- make.tree.classe.td(pars.cl, k=3, max.t1=30, max.t2=35, x0=1, single.lineage=TRUE)
tb
phy <- table2tree(tb)
phy$epoch1.ntip
Ntip(phy)- phy$epoch1.ntip
plot(phy)

phy$t.epoch1
phy$t.total
phy$epoch1.extant.tree.depth
phy$epoch12.tree.depth
phy$epoch1.ntip
Ntip(phy)
#phy$tip.state
#length(phy$tip.state)
