library(divextra)
library(plyr)
source('R/asr.R')


#------- Toy dataset, parameters given
data(phy5)
cols <- plyr::mapvalues(phy5$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
plot(phy5, label.offset = 2, cex=1, no.margin = F)
tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
nodelabels()
tiplabels()
edgelabels()
axisPhylo()

lik.td <-make.classe.td(phy5, phy5$tip.state_3, k=3, n.epoch=2, control=list(backend = "gslode"), strict=T)
par.td <- c(4, c(1:54)/100)
names(par.td) <- diversitree::argnames(lik.td)
lik.est <- lik.td(par.td, intermediates=T)
print(lik.est)

st <-  asr.marginal.classe(lik.td, par.td)
nodelabels(thermo=t(st), piecol=c('black', "blue","red"), cex=2, adj=.5)


#------- Toy dataset, parameters estimated
file_path <- system.file("extdata", "geosse_tb_tree.rds", package = "divextra")
phy <- readRDS(file_path)
plot(phy)

file_path <- system.file("extdata", "geosse3_td_test.yml", package = "divextra")
par.categories.td <- read_yaml_pars_td(file_path)
# display
print(par.categories.td)

lik.td <-make.classe.td(phy, phy$tip.state, k=3, n.epoch=2, control=list(backend = "gslode"), strict=T)
formula.td <- make_constraints_sse_td(par.categories.td)
lik.const.td <- constrain(lik.td, formulae = formula.td)
starting.point <- init.pars.classe_td(lik.const.td, phy, k=3, n.epoch=2, eps=0.5)
print(starting.point)

mle.td <- find.mle(lik.const.td, starting.point, condition.surv=TRUE, keep.func=F)

st <- asr.marginal.classe(lik.td, mle.td$par.full)

cols <- plyr::mapvalues(phy$tip.state, from = c(1:3), to=c('orange', "green","red" ))
plot(phy, show.tip.label = F, cex=1, no.margin = F)
tiplabels(pch = 15, col = cols, cex = .5, adj = 1)
nodelabels(thermo=t(st), piecol=c('orange', "green","red" ), cex=.3, adj=0)


#------------------------------
#---- ASR
data(phy5)
cols <- plyr::mapvalues(phy5$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
plot(phy5, label.offset = 2, cex=1, no.margin = F)
tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
nodelabels()
tiplabels()
edgelabels()
axisPhylo()

lik.td <-make.classe.td(phy5, phy5$tip.state_3, k=3, n.epoch=2, control=list(backend = "gslode"), strict=T)
par.td <- c(4, c(1:54)/100)
names(par.td) <- diversitree::argnames(lik.td)
lik.est <- lik.td(par.td, intermediates=T)
print(lik.est)

st <-  asr.marginal.classe(lik.td, par.td)
# ASR reconstruction at t=1 of the branch 3 (defined by the terminal node 9)
# Note: relative.t.branch starts from present time and rquals to 0 at the branch's start
st9 <-asr.marginal.classe_branch(lik.td, par.td, node.id=9, relative.t.branch= 1 )

print(st)
print(st9)

#------------------------------
#---- ASR
# make.asr.marginal.classe(lik.td)(par.td, nodes=4)
st <- make.asr.marginal.classe(lik.td)(par.td)
st
st <- asr.marginal.classe(lik.td, par.td)
# cols <- mapvalues(phy$tip.state_3, from = c(1:3), to=c('black', "blue","red" ))
# plot(phy, label.offset = 2, cex=1, no.margin = F)
# tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
# nodelabels(thermo=t(st), piecol=c('black', "blue","red"), cex=2, adj=.5)

#-----




#---- ASR along branch
node.plot <- 9
node.id <- node.plot - 5
node.id
#relative.t.branch
source('R/asr.R')
phy$edge.length # 5.564629 # depth: 1.502053
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= 0)
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= 5.564)
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= phy$edge.length[1])
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= "BASE")
st
asr.marginal.classe_branch.multiple(lik.td, par.td, node.id=9, Nbins=10, include.ends=TRUE)

asr.marginal.classe_branch(lik.td, par.td, node.id=9, relative.t.branch= 5.9)
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= 5)
asr.marginal.classe_branch(lik.td, par.td, node.id=9, relative.t.branch= "BASE")
make.asr.marginal.classe_branch(lik.td)(par.td, nodes=node.id, relative.t.branch= "BASE")
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




