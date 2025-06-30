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

#tol=1e-16, eps=1e-4
lik.td <-make.classe.td(phy, phy$tip.state_3, k=3, n.epoch=2,
                        control=list(backend = "gslode"), strict=T)
par.td <- c(4, c(1:54)/100)
names(par.td) <- diversitree::argnames(lik.td)
lik.est <- lik.td(par.td, intermediates=T)
#print(lik.est)

#---- ASR
# make.asr.marginal.classe(lik.td)(par.td, nodes=4)
st <- make.asr.marginal.classe(lik.td)(par.td, root=ROOT.GIVEN, root.p=c(1,0,0))
st
# cols <- mapvalues(phy$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
# plot(phy, label.offset = 2, cex=1, no.margin = F)
# tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
# nodelabels(thermo=t(st), piecol=c('black', "blue","red"), cex=2, adj=.5)



#---- ASR along branch
node.plot <- 7
node.id <- node.plot - 5
node.id
#relative.t.branch
source('R/asr.R')
phy$edge.length # 5.564629 # depth: 1.502053
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= 0, root=ROOT.GIVEN, root.p=c(1,0,0))
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= 5.564)
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= phy$edge.length[1])
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= "BASE")
st

# [,1]
# [1,] 0.3646077
# [2,] 0.3308331
# [3,] 0.3045592
#--- ASR along branch
Nbins <- 10
lik <- lik.td
eps <- 0.01
e <- environment(lik)
eb <- environment(e$all_branches)
cache <- eb$cache

node.plot <- 9 # upstream node
br.len <- cache$len[node.plot]
node.id.diversitree <- node.plot - Ntip(phy)

time.bins <- seq(0+eps, br.len-eps, length.out=Nbins)

ans.t <- c()
for (t in time.bins){
  print(t)
  x <- make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id.diversitree,
                                               relative.t.branch= t, root=ROOT.GIVEN, root.p = c(1,0,0))
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

#---- ASR along branch. Markdown

st <- make.asr.marginal.classe(lik.td)(mle.td$par.full)
#st
#round(st, 2)
round(st[,1:10], 2)

cols <- mapvalues(phy$tip.state, from = c(1:3), to=c('orange', "green","red" ))
plot(phy, show.tip.label = F, cex=1, no.margin = F)
nodelabels(cex=.3)
tiplabels(pch = 15, col = cols, cex = .5, adj = 1)
nodelabels(thermo=t(st), piecol=c('orange', "green","red" ), cex=.3, adj=0)
axisPhylo()

Nbins <- 10
lik <- lik.td
eps <- 0.01
e <- environment(lik)
eb <- environment(e$all_branches)
cache <- eb$cache

node.plot <- 801 # upstream node
br.len <- cache$len[node.plot]
node.id.diversitree <- node.plot - Ntip(phy)

time.bins <- seq(0+eps, br.len-eps, length.out=Nbins)

ans.t <- c()
for (t in time.bins){
  print(t)
  x <- make.asr.marginal.classe_branch(lik.td)(mle.td$par.full, nodes=node.id.diversitree,
                                               relative.t.branch= t)
  ans.t <- cbind(ans.t, x)
}

st[,node.id.diversitree]
st[,729-Ntip(phy)]
ans.t
plot(time.bins, ans.t[3,], type='l')
