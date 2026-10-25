library(divextra)

# Rows are oldest to most recent ----
pars.ge <- rbind(
  c(.1,.1,.1,.05,.05,.1,.1),
  c(.3,.3,.3,.05,.05,.1,.1))
pars.cl <- rbind(pars.geosse2classe(pars.ge[1,]),
                 pars.geosse2classe(pars.ge[2,]))

# Switch when 100 living lineages are reached, or after 60 time units.
# Stop the second epoch at 500 living lineages, or after another 10 units.
set.seed(10)
tb <- make.tree.classe.episodic(
  pars.cl, k=3, x0=1,
  ntips=c(100,500), duration=c(60,10))

attr(tb,"epochs")
attr(tb,"completed")
attr(tb,"stop.reason")

# Boundary age for a two-epoch likelihood; not a difference of crown depths.
boundary.age <- attr(tb,"epochs")$end.age[1]

# Histories include extinct lineages. A completely unsplit history may not
# convert to a branching phylogeny; inspect completed and the table first.
if(sum(tb$split)>0) {
  mapped <- table2simmap(tb,state.labels=c("A","B","AB"))
  phytools::plotSimmap(mapped,colors=c(A="#0072B2",B="#E69F00",AB="#009E73"))
}
