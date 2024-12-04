library(yaml)
library(diversitree)
library(divextra)
source('R/utils-yaml.R')
source('R/utils-sse.R')

# make.musse.td
# make.bisse.td
# make.classe

# library(diversitree)

## GeoSSE equivalence
## Same tree simulated in ?make.geosse
# pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
# names(pars) <- diversitree:::default.argnames.geosse()
# set.seed(5)
# phy <- tree.geosse(pars, max.taxa=5, x0=0)
# saveRDS(phy, 'my-test/data/phy.rds')
# write.tree(phy, 'my-test/data/phy.tree')

## GeoSSE equivalence
## Same tree simulated in ?make.geosse
# pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
# names(pars) <- diversitree:::default.argnames.geosse()
# set.seed(5)
# phy <- tree.geosse(pars, max.taxa=50, x0=0)
# saveRDS(phy, 'my-test/data/phy-50tips.rds')
# write.tree(phy, 'my-test/data/phy.tree')

#------
phy <- readRDS('my-test/data/phy.rds')
plot(phy)

lik.g <- make.geosse(phy, phy$tip.state)
pars.g <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
names(pars.g) <- argnames(lik.g)


lik.c <- make.classe(phy, phy$tip.state+1, 3)
pars.c <- starting.point.classe(phy, 3)

# pars.c <- 0 * starting.point.classe(phy, 3)
# pars.c['lambda222'] <- pars.c['lambda112'] <- pars.g['sA']
# pars.c['lambda333'] <- pars.c['lambda113'] <- pars.g['sB']
# pars.c['lambda123'] <-  pars.g['sAB']
# pars.c['mu2'] <- pars.c['q13'] <- pars.g['xA']
# pars.c['mu3'] <- pars.c['q12'] <- pars.g['xB']
# pars.c['q21'] <- pars.g['dA']
# pars.c['q31'] <- pars.g['dB']

lik.g(pars.g)   # xx
lik.c(pars.c)   # -8.367568

#----------
## For comparison, make a plain MuSSE likelihood function
lik.m <- make.musse(phy, phy$tip.state+1, 3)

## Create the time-dependent likelihood function.  The final argument
## here is the number of 'epochs' that are allowed.  Two epochs is one
## switch point.
pars <- c(.1,  .15,  .2,  # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0,         # q12, q13
          .05, .05,       # q21, q23
          0,   .05)       # q31, q32
lik.t <- make.musse.td(phy, phy$tip.state+1, 3, 2)
argnames(lik.t)
pars.t <- c(3, pars, pars)
names(pars.t) <- argnames(lik.t)

lik.t(pars.t)
lik.t2 <- constrain(lik.t, t.1 ~ 3)
## And fit the MuSSE/td model
fit.t <- find.mle(lik.t2, pars.t[argnames(lik.t2)],control=list(maxit=20000))

#--- Classe td

lik.ctd <- make.classe.td(phy, phy$tip.state+1, k=3, n.epoch=2, control=list(backend = "gslode"))

phy$edge.length
plot(phy)
edgelabels(phy$edge.length)

argnames(lik.ctd)
par.ctd <- c(0.2, pars.c, pars.c)
names(par.ctd) <- argnames(lik.ctd)
par.ctd
#par.ctd['lambda222.2'] <- 0.4
#par.ctd['lambda222.1'] <- 0.1
lik.ctd(par.ctd)
#-8.367568



argnames(lik.ctd)
par.ctd <- c(0.2, pars.c, pars.c)
names(par.ctd) <- argnames(lik.ctd)
par.ctd
par.ctd['lambda222.2'] <- 0.4
par.ctd['lambda222.1'] <- 0.1
lik.ctd(par.ctd)
#-8.367568

#----------------------------------------
phy <- read.tree(text = '(sp1:1.014858647,(((sp4:0.2614280769,sp5:0.2614280769)nd4:0.007869378189,sp3:0.2692974551)nd3:0.01878075832,sp2:0.2880782134)nd2:0.7267804339)nd1;')
states <- species_counts <- c(sp1 = 2, sp2 = 2, sp3 = 1, sp4 = 0, sp5 = 1)

trace("make.classe", tracer = quote(cat("function_name called\n")),  where = "package:diversitree")

debug(diversitree::make.classe)
debug(lik.c)

# classe model
lik.c <- diversitree::make.classe(phy, states+1, 3)
pars.c <- starting.point.classe(phy, 3)
lik.c(pars.c)   # -8.367568




# classe.td where two time bins have the same pars (=classe)
lik.ctd <- make.classe.td(phy, states + 1, k=3, n.epoch=2, control=list(backend = "gslode"))
par.ctd <- c(0.2, pars.c, pars.c)
names(par.ctd) <- argnames(lik.ctd)
par.ctd
lik.ctd(par.ctd) # -8.367568

# classe.td where two time bins have different pars
lik.ctd2 <- make.classe.td(phy, states + 1, k=3, n.epoch=2, control=list(backend = "gslode"))
par.ctd2 <- c(0.2, pars.c, pars.c)
names(par.ctd2) <- argnames(lik.ctd)
par.ctd2['lambda222.2'] <- 0.4
par.ctd2['lambda222.1'] <- 0.1
lik.ctd(par.ctd2) # -8.268788


#---------------------------------------- Prune node
phy <- read.tree(text = '(sp1:1.014858647,(((sp4:0.2614280769,sp5:0.2614280769)nd4:0.007869378189,sp3:0.2692974551)nd3:0.01878075832,sp2:0.2880782134)nd2:0.7267804339)nd1;')


library(ape)

plot(phy)
nodelabels()

node_to_prune <- 9

# Get all descendants of the specified node
descendants <- phangorn::Descendants(phy, node = node_to_prune, type = "all")

# Prune the descendants
pruned_tree <- drop.tip(phy, phy$tip.label[descendants], trim.internal = FALSE)

# Plot the pruned tree
plot(pruned_tree, main = "Tree After Pruning Descendants")
nodelabels()

# classe model
states.pruned <- species_counts <- c(sp1 = 3, sp2 = 2, sp3 = 1, nd4=1)
lik.c <- make.classe.td(pruned_tree, states.pruned, 3, strict=F)
pars.c <- starting.point.classe(phy, 3)
lik.c(pars.c)   # -8.367568



st.m <- asr.marginal(lik, pars[5:6])
