library(divextra)

if (file.exists("R/sim-on-tree-classe-td.R")) {
  if (!requireNamespace("pkgload", quietly = TRUE))
    stop("Install pkgload or install the current DivExtra checkout first")
  pkgload::load_all(".", export_all = FALSE)
} else {
  library(divextra)
}

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
tb <- make.tree.classe.td(pars.cl, k=3, max.t1=10, max.t2=15, x0=2, single.lineage=TRUE)
print(tb)
phy <- table2tree(tb)
plot(phy)

simmap <- table2simmap(
  tb,
  state.labels = c("A", "B", "AB")
)


colors <- c(
  A  = "green",
  B  = "#0072B2",
  AB = "red"
)

phytools::plotSimmap(
  simmap,
  colors = colors,
  ftype = "i",
  split.vertical = TRUE
)

#----- make html

library(divextra)

## GeoSSE parameter order:
## "sA" "sB" "sAB" "xA" "xB" "dA" "dB"

pars.ge <- matrix(
  c(
    0.1, 0.1, 0.1, 0, 0, 0.1, 0.1,
    0.3, 0.3, 0.3, 0, 0, 0.1, 0.1
  ),
  nrow = 2,
  ncol = 7,
  byrow = TRUE
)

colnames(pars.ge) <- diversitree:::default.argnames.geosse()

# Convert each epoch from GeoSSE to ClaSSE parameters
pars.cl <- t(apply(pars.ge, 1, pars.geosse2classe))

# Add epoch suffixes required by create_parameter_html()
parameter_values <- unlist(
  lapply(seq_len(nrow(pars.cl)), function(epoch) {
    setNames(
      pars.cl[epoch, ],
      paste0(colnames(pars.cl), ".", epoch)
    )
  })
)

# Zero-valued parameters are implicit and should be omitted
parameter_values <- parameter_values[parameter_values != 0]

# Construct the object expected by create_parameter_html()
par.categories.td <- list(
  Nstates    = 3L,
  states     = c("A", "B", "A.B"),
  n.epoch    = nrow(pars.cl),
  epoch.times = 10,  # replace with your actual epoch boundary
  pars       = split(
    names(parameter_values),
    as.character(parameter_values)
  )
)

# Ensure the output directory exists
##dir.create("sim-on-tree", recursive = TRUE, showWarnings = FALSE)

create_parameter_html(
  par.categories.td = par.categories.td,
  model_name = "GeoSSE Test",
  output_file = "sim-on-tree/geosse-test.html"
)


#----------






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
ape::Ntip(phy) # total N of tips

# The model changes regimes at (if counting from the tips for the infference):
phy$t.regime.change
