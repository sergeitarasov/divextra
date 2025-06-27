library(diversitree)
library(plyr)

# pars <- c(.1, .2, .03, .06, .01, .02)
# phy <- trees(pars, "bisse", max.taxa=5, max.t=Inf, x0=0)[[1]]
# phy$tip.state <- c(1,1, 0,0,1)
# names(phy$tip.state) <- phy$tip.label
# saveRDS(phy, 'asr/tree.rds')
phy <- readRDS('asr/tree.rds')


cols <- mapvalues(phy$tip.state, from = c(0, 1), to=c("grey", 'black'))
plot(phy, label.offset = 2, cex=1, no.margin = F)
tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
nodelabels()
tiplabels()
edgelabels()
axisPhylo()

pars <- c(.1, .2, .03, .06, .01, .02)

lik <- make.bisse(phy, phy$tip.state)
st <- asr.marginal(lik, pars, root=ROOT.FLAT)
plot(phy, label.offset = 2, cex=1, no.margin = F)
tiplabels(pch = 15, col = cols, cex = 2, adj = 1.5)
nodelabels(thermo=t(st), piecol=c(0,1), cex=3)

#----------------


#-----
diversitree:::make.asr.marginal.bisse(lik, pars)

asr.bisse <- diversitree:::make.asr.marginal.bisse(lik)
asr.bisse(pars)

diversitree:::make.asr.marginal.bisse(lik)(pars)

#-------------------




