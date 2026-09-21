#library(devtools)
#install_github("sergeitarasov/divextra")
library(ape)
library(phytools)
library(readxl)
library(geiger)
library(plyr)
library(dplyr)
library(phangorn)
library(tidyverse)
library(jsonlite)
library(rstudioapi)
library(fs)
library(divextra)


# this returns the full path to the active document
# active_path <- getActiveDocumentContext()$path
# active_dir <- dirname(active_path)
# setwd(active_dir)

#-----------------------------------
# PARAMETERS
# the params in {{}} are generated atomatically
#-----------------------------------
# PARAM_GR=1
# RUN =1
PARAM_GR= 6
RUN = 3
BASE_NAME="AfrSunken16"
YAML_BACKBONE = "sim-on-tree/yml/backbone-Afr-Sunken-16.yml"
YAML_GROUPS = "sim-on-tree/yml/groups-Afr-Sunken-16.yml"
EPSs <- seq(0.1, 0.8, length.out=5)

# #----------------------------------
EPS=EPSs[RUN]
OUTPUT_BASE <- paste0(BASE_NAME, '_gr-', PARAM_GR)
print(OUTPUT_BASE)
# print(RUN)
#
# #-----------------------------------
# # Read Data
# #-----------------------------------
phy <- readRDS("sim-on-tree/phy-271.rds")
geo_regions_abc <- readRDS("sim-on-tree/geo_regions-271.rds")
#
# states: [A, E, M, U, S, R,  ...]
geo_regions <- mapvalues(geo_regions_abc, from = c("Afr", "OP", "Mada", "Aus", "Maur"), to=c(1:4, 6) ) %>% as.numeric()
names(geo_regions) <- names(geo_regions_abc)
#
# #-----------------------------------
# # sampling fraction
# #-----------------------------------
# # states:                [A,       E,     M,       U,      S,      R,  ...]
single_transitions.f = c(0.039,   0.027,  0.446, 0.008,   1e-06,   0.833)
dual_transitions.f <- rep(1e-06, 10)
sampling.f <- c(single_transitions.f, dual_transitions.f)
names(sampling.f) <- 1:16
sampling.f
#
# #-----------------------------------
# # Read Backbone Params
# #-----------------------------------
par.categories.td <- read_yaml_pars_td(YAML_BACKBONE)
print(par.categories.td)
#
# #-----------------------------------
# # Read Params Groups
# #----------------------------------
par_groups_yaml <- read_par_groups_yaml(YAML_GROUPS)
par_groups_yaml
#
# #-----------------------------------
# #  Create Focal Params
# #----------------------------------
focal_group <- par_groups_yaml[[PARAM_GR]]
focal_group
new_pars <- regroup_parameters(par.categories.td, focal_group)

# #-----------------------------------
# # ML setup
# #-----------------------------------
mle.td <- NULL

lik.td <-make.classe.td(phy, geo_regions, k=16, n.epoch=2, control=list(backend ="gslode"), strict=F, sampling.f=sampling.f)
formula.td <- make_constraints_sse_td(new_pars)
lik.const.td <- constrain(lik.td, formulae = formula.td)
starting.point <- init.pars.classe_td(lik.const.td, phy, k=16, n.epoch=2, eps= EPS)
print(starting.point)



#-----------------------------------
# ASR
#-----------------------------------

mle.td <- readRDS(file=path("sim-on-tree", paste0(OUTPUT_BASE, "_r-", RUN, ".rds")))
mle.td$AIC

st <- asr.marginal.classe(lik.td, mle.td$par.full, root=ROOT.GIVEN, root.p=c(1, rep(0,15)))
rownames(st) <- c('A', 'E', 'M', 'U', 'S', 'R',   'A.E', 'A.M', 'A.U', 'A.S', 'A.R',   'E.M', 'E.U',  'E.S', 'E.R',    'S.R')
round(st[,1:10], 2)

# states: [A, E, M, U, S, R,  ...]

# Plot
cols <- mapvalues(geo_regions_abc, from = c("Aus", "Mada", "OP", "Afr", "Maur"), to=c("yellow","blue", "green","red", "purple") )
plot(phy, 'p', label.offset = .008, cex=0.2, no.margin = F)
tiplabels(pch = 15, col = cols, cex = .3, adj = .503)
nodelabels(cex = .3, height = .05, width = .05, frame='none', col ='red')
axisPhylo()


#-----------


# states: [A, E, M, U, !S, R,   A.E, A.M, A.U, A.S, A.R,   E.M, E.U,    S.R]
colsDB <- c(
  # Single regions
  "A"   = "#D73027",  # Strong red
  "E"   = "#A6D854",  # Rich green
  "M"   = "#74ADD1",  # Vivid pink
  "U"   = "grey",  # Bright yellow
  "S"   = "#E78AC3",  # Deep blue
  "R"   = "#FFD92F",  # Light teal-blue (distinct from R)

  # Combinations involving A
  "A.E" = "#6DAF27",  # Brick red
  "A.M" = "#B2182B",  # Coral
  "A.U" = "#F4A582",  # Soft orange
  "A.S" = "#F46D43",  # Salmon-orange
  "A.R" = "white",  # Soft blue (Africa + Mauritius)

  # Eurasian combos
  "E.M" = "#66C2A5",  # Muted teal
  "E.U" = "white",  # Lime green
  'E.S' = "#66C2A1",
  'E.R' = "#66C2A9",
  "S.R" = "black"   # Muted periwinkle-blue (for balance between S and R)
)

colsDB
# plotTree(phy, type='fan', fsize=0.1,ftype="i")
plotTree(phy,fsize=0.1, ftype="i", lwd=1, offset=1)
nodelabels(node=1:phy$Nnode+Ntip(phy), pie=t(st), piecol=colsDB, cex=0.3)


# tiplabels(pie=to.matrix(geo_regions,sort(unique(geo_regions))), piecol=colsDB[-5], cex=0.2, adj = 2)
# plot()
# legend("bottomleft",
#        legend = names(colsDB),
#        fill = colsDB,
#        border = NA,
#        bty = "n",     # no box
#        cex = 1.5,     # text size
#        pt.cex = 1.2)  # square size



# pdf("runs/asr/asr-pres-zoom.pdf", width = 20, height = 12)  # Change size as needed
plotTree(phy, fsize=0.3, ftype="i", lwd=2, offset=.2)
nodelabels(node=1:phy$Nnode+Ntip(phy), pie=t(st), piecol=colsDB, cex=0.11)
tiplabels(pie=to.matrix(geo_regions,sort(unique(geo_regions))), piecol=colsDB[-5], cex=0.05, adj = 0.5)
# dev.off()

#---------------------------------
#-----------  Extract Neosisyphus
extract.clade.simm<-function(tree,node){
  x<-getDescendants(tree,node)
  x<-x[x<=length(tree$tip.label)]
  drop.tip(tree,tree$tip.label[-x])
}
# ASR plot
x<-getDescendants(phy, 495)
#Descendants(phy, 495, type = c( "children"))
x<-x[x >length(phy$tip.label)]
length(x)
x <- c(495, x)
x <- x[order(x)]
#Descendants(x, node, type = c("tips", "children", "all"))
phySis <- extract.clade.simm(phy, 495)
stSis <- st[,x-Ntip(phy)]

geo_regionsSis <- geo_regions[match(phySis$tip.label, names(geo_regions))]

#pdf("runs/asr/asr-Nesosisiphus.pdf", width = 10, height = 6)  # Change size as needed
plotTree(phySis, fsize=0.7, ftype="i", lwd=3, offset=.6)
nodelabels(node=1:phySis$Nnode+Ntip(phySis), pie=t(stSis), piecol=colsDB, cex=1)
tiplabels(pie=to.matrix(geo_regionsSis,sort(unique(geo_regionsSis))), piecol=colsDB[-c(3,4,5)], cex=0.4, adj = 0.5)

