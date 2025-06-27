library(divextra)
library(plyr)
source('R/asr.R')

# Read in tree and Plot
phy <- readRDS('asr/tree.rds')
cols <- mapvalues(phy$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
plot(phy, label.offset = 2, cex=1, no.margin = F)
tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
nodelabels()
nodelabels(node=c(7))
tiplabels()
edgelabels()
axisPhylo()

lik.td <-make.classe.td(phy, phy$tip.state_3, k=3, n.epoch=2, control=list(backend = "gslode", tol=1e-16, eps=1e-4), strict=T)
par.td <- c(4, c(1:54)/100)
names(par.td) <- diversitree::argnames(lik.td)
lik.est <- lik.td(par.td, root=ROOT.FLAT, intermediates=T)
print(lik.est)

#---- ASR
make.asr.marginal.classe(lik.td)(par.td, nodes=4)
st <- make.asr.marginal.classe(lik.td)(par.td)
st
cols <- mapvalues(phy$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
plot(phy, label.offset = 2, cex=1, no.margin = F)
tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
nodelabels(thermo=t(st), piecol=c('black', "blue","red"), cex=2, adj=.5)

getMRCA(phy, c(1,2))

#---- ASR along branch
node.plot <- 9
node.id <- node.plot - 5
node.id
#relative.t.branch
source('R/asr.R')
phy$edge.length # 5.564629
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= 0.001)
# 0.1822480

#--- ASR along branch
Nbins <- 10
lik <- lik.td
eps <- 0.01
e <- environment(lik)
eb <- environment(e$all_branches)
cache <- eb$cache

node.plot <- 8 # upstream node
br.len <- cache$len[node.plot]
node.id.diversitree <- node.plot - Ntip(phy)

time.bins <- seq(0+eps, br.len-eps, length.out=Nbins)

ans.t <- c()
for (t in time.bins){
  print(t)
  x <- make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id.diversitree, relative.t.branch= t)
  ans.t <- cbind(ans.t, x)
}

st
ans.t
plot(time.bins, ans.t[3,], type='l')

# nodes + cache$n.tip
#
# children <- cache$children
# parent <- cache$parent
# len <- cache$len
# depth <- cache$depth
# root.idx <- cache$root
# anc <- cache$ancestors




