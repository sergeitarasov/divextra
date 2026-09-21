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
RUN = 1
BASE_NAME="AfrSunken"
YAML_BACKBONE = "yml/backbone-Sunken.yml"
YAML_GROUPS = "yml/groups-Sunken.yml"
EPSs <- seq(0.1, 0.8, length.out=5)

#----------------------------------
EPS=EPSs[RUN]
OUTPUT_BASE <- paste0(BASE_NAME, '_gr-', PARAM_GR)
print(OUTPUT_BASE)
print(RUN)

#-----------------------------------
# Read Data
#-----------------------------------
phy <- readRDS("data/phy-271.rds")
geo_regions_abc <- readRDS("data/geo_regions-271.rds")

# states: [A, E, M, U, S, R,  ...]
geo_regions <- mapvalues(geo_regions_abc, from = c("Afr", "OP", "Mada", "Aus", "Maur"), to=c(1:4, 6) ) %>% as.numeric()
names(geo_regions) <- names(geo_regions_abc)

#-----------------------------------
# sampling fraction
#-----------------------------------
# states:                [A,       E,     M,       U,      S,      R,  ...]
single_transitions.f = c(0.039,   0.027,  0.446, 0.008,   1e-06,   0.833)
dual_transitions.f <- rep(1e-06, 8)
sampling.f <- c(single_transitions.f, dual_transitions.f)
names(sampling.f) <- 1:14
sampling.f

#-----------------------------------
# Read Backbone Params
#-----------------------------------
par.categories.td <- read_yaml_pars_td(YAML_BACKBONE)
print(par.categories.td)

#-----------------------------------
# Read Params Groups
#----------------------------------
par_groups_yaml <- read_par_groups_yaml(YAML_GROUPS)
par_groups_yaml

#-----------------------------------
#  Create Focal Params
#----------------------------------
focal_group <- par_groups_yaml[[PARAM_GR]]
focal_group
#focal_group_name <- names(par_groups_yaml[1])
new_pars <- regroup_parameters(par.categories.td, focal_group)
#new_pars
# create_parameter_html(
#   par.categories.td = new_pars,
#   model_name = OUTPUT_BASE,
#   output_file =  path("../", "html", paste0(OUTPUT_BASE, ".html"))
# )


#-----------------------------------
# ML setup
#-----------------------------------
mle.td <- NULL

lik.td <-make.classe.td(phy, geo_regions, k=14, n.epoch=2, control=list(backend ="gslode"), strict=F, sampling.f=sampling.f)
formula.td <- make_constraints_sse_td(new_pars)
lik.const.td <- constrain(lik.td, formulae = formula.td)
starting.point <- init.pars.classe_td(lik.const.td, phy, k=14, n.epoch=2, eps= EPS)
print(starting.point)

#mle.td <- find.mle(lik.const.td, starting.point, condition.surv=TRUE, keep.func=F, root=ROOT.OBS, control=list(maxit=20000))
mle.td <- tryCatch(
  find.mle(lik.const.td, starting.point, condition.surv=TRUE, keep.func=F, root=ROOT.GIVEN, root.p=c(1, rep(0,13)), control=list(maxit=20000)),
  error = function(e) NULL
)

# Simple check and continue
if (is.null(mle.td)) {
  message(sprintf("Failed at EPS=%s", EPS))
  next
}


#-----------------------------------
# ASR
#-----------------------------------

mle.td <- readRDS(file=path("runs", "rds", paste0(OUTPUT_BASE, "_r-", RUN, ".rds")))
mle.td$AIC

st <- asr.marginal.classe(lik.td, mle.td$par.full, root=ROOT.GIVEN, root.p=c(1, rep(0,13)))
rownames(st) <- c('A', 'E', 'M', 'U', 'S', 'R',   'A.E', 'A.M', 'A.U', 'A.S', 'A.R',   'E.M', 'E.U',    'S.R')
round(st[,1:10], 2)

# states: [A, E, M, U, S, R,  ...]
#unique(geo_regions_abc)
#table(geo_regions_abc)
# geo_regions <- mapvalues(geo_regions_abc, from = c("Afr", "OP", "Mada", "Aus", "Maur"), to=c(1:4, 6) ) %>% as.numeric()
# #unique(geo_regions)
# names(geo_regions) <- names(geo_regions_abc)

# Plot
cols <- mapvalues(geo_regions_abc, from = c("Aus", "Mada", "OP", "Afr", "Maur"), to=c("yellow","blue", "green","red", "purple") )
plot(phy, 'p', label.offset = .008, cex=0.2, no.margin = F)
tiplabels(pch = 15, col = cols, cex = .3, adj = .503)
nodelabels(cex = .3, height = .05, width = .05, frame='none', col ='red')
axisPhylo()


br <- asr.marginal.classe_branch.multiple(lik.td, mle.td$par.full, node.id=513, eps = 0.1, Nbins=10, include.ends=TRUE)
round(st[,1:10], 2)
round(t(st[,1:10]), 2)

# states: [A, E, M, U, S, R,   A.E, A.M, A.U, A.S, A.R,   E.M, E.U,    S.R]
plot(br$t + br$depth, br$asr[10,], type='l', xlim=rev(range(br$t)))
plot(br$t + br$depth, br$asr[5,], type='l', xlim=rev(range(br$t)))

513-Ntip(phy)
495-Ntip(phy)
st[,242]
st[,224]

st[,513-Ntip(phy)]
phy$edge.length[,272]
node.depth.edgelength(phy)[272]
max(node.depth.edgelength(phy))

node.depth.edgelength(phy)[495] - 96.98
node.depth.edgelength(phy)[489] - 96.98
#-----------

colsDB <- c(
  # Single regions
  "A"   = "#D73027",  # Strong red
  "E"   = "#A6D854",  # Rich green
  "M"   = "#74ADD1",  # Vivid pink
  "U"   = "grey",  # Bright yellow
  "S"   = "#E78AC3",  # Deep blue
  "R"   = "#FFD92F",  # Light teal-blue (distinct from R)

  # Combinations involving A
  "A.E" = "#B2182B",  # Brick red
  "A.M" = "#EF8A62",  # Coral
  "A.U" = "#F4A582",  # Soft orange
  "A.S" = "#FDAE6B",  # Salmon-orange
  "A.R" = "#A6CEE3",  # Soft blue (Africa + Mauritius)

  # Eurasian combos
  "E.M" = "#66C2A5",  # Muted teal
  "E.U" = "#A6D854",  # Lime green
  "S.R" = "#80B1D3"   # Muted periwinkle-blue (for balance between S and R)
)

# states: [A, E, M, U, !S, R,   A.E, A.M, A.U, A.S, A.R,   E.M, E.U,    S.R]
apply(st, 1, sum)
apply(st, 1, max)


"#FDAE6B"
"#F46D43"
"#F4A6C0"

"#C967A3"
"#DB91B2"

original <- "#E78AC3"  # existing color (pink)
similar  <- "#D47DBB"  # distinct enough to sit nearby

"base"     = "#A6D854",
"brighter" = "#C5E887",
"lighter"  = "#DFF3B4",
"darker"   = "#6DAF27",
"deeper"   = "#7CB342",
"muted"    = "#A9C97E",
"vibrant"  = "#9CD700"
"#66C2A5"

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
  "S.R" = "black"   # Muted periwinkle-blue (for balance between S and R)
)

colsDB
#colsDB<-setNames(palette()[1:length(unique(geo_regions))],sort(unique(geo_regions)))

plotTree(phy, type='fan', fsize=0.1,ftype="i")

pdf("runs/asr/asr-pres.pdf", width = 10, height = 6)  # Change size as needed
plotTree(phy,fsize=0.1, ftype="i", lwd=1, offset=1)
nodelabels(node=1:phy$Nnode+Ntip(phy), pie=t(st), piecol=colsDB, cex=0.3)
#nodelabels(cex = .3, height = .05, width = .05, frame='none', col ='grey')
tiplabels(pie=to.matrix(geo_regions,sort(unique(geo_regions))), piecol=colsDB[-5], cex=0.2, adj = 2)
#axisPhylo()
# --- Close device ---
# Add legend (adjust placement as needed)
plot()
legend("bottomleft",
       legend = names(colsDB),
       fill = colsDB,
       border = NA,
       bty = "n",     # no box
       cex = 1.5,     # text size
       pt.cex = 1.2)  # square size
dev.off()


pdf("runs/asr/asr-pres-zoom.pdf", width = 20, height = 12)  # Change size as needed
plotTree(phy, fsize=0.3, ftype="i", lwd=2, offset=.2)
nodelabels(node=1:phy$Nnode+Ntip(phy), pie=t(st), piecol=colsDB, cex=0.11)
#nodelabels(cex = .3, height = .05, width = .05, frame='none', col ='grey')
tiplabels(pie=to.matrix(geo_regions,sort(unique(geo_regions))), piecol=colsDB[-5], cex=0.05, adj = 0.5)
#axisPhylo()
# --- Close device ---
dev.off()

#---------------------------------
#----------- Neosisyphus
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

pdf("runs/asr/asr-Nesosisiphus.pdf", width = 10, height = 6)  # Change size as needed
plotTree(phySis, fsize=0.7, ftype="i", lwd=3, offset=.6)
nodelabels(node=1:phySis$Nnode+Ntip(phySis), pie=t(stSis), piecol=colsDB, cex=1)
#nodelabels(cex = .3, height = .05, width = .05, frame='none', col ='grey')
tiplabels(pie=to.matrix(geo_regionsSis,sort(unique(geo_regionsSis))), piecol=colsDB[-c(3,4,5)], cex=0.4, adj = 0.5)
#axisPhylo()
# --- Close device ---
dev.off()

#-------------

br <- asr.marginal.classe_branch.multiple(lik.td, mle.td$par.full, node.id=513, eps = 0.1, Nbins=100, include.ends=TRUE)

# states: [A, E, M, U, S, R,   A.E, A.M, A.U, A.S, A.R,   E.M, E.U,    S.R]

# A
plot(br$t + br$depth, br$asr[1,], type='l', xlim=rev(range(br$t)))
# A.S
plot(br$t + br$depth, br$asr[10,], type='l', xlim=rev(range(br$t)))
# S
plot(br$t + br$depth, br$asr[5,], type='l', xlim=rev(range(br$t)))
# S.R
plot(br$t + br$depth, br$asr[14,], type='l', xlim=rev(range(br$t)))
# R
plot(br$t + br$depth, br$asr[6,], type='l', xlim=rev(range(br$t)))

#--------
# A
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[1,], type='l', xlim=rev(range(brNesovin$t)))
lines(br$t + br$depth, br$asr[1,], type='l')

# A.S
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[10,], type='l', xlim=rev(range(brNesovin$t)))
lines(br$t + br$depth, br$asr[10,], type='l')

#S
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[5,], type='l', xlim=rev(range(brNesovin$t)))
lines(br$t + br$depth, br$asr[5,], type='l')

# S.R
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[14,], type='l', xlim=rev(range(brNesovin$t)))
lines(br$t + br$depth, br$asr[14,], type='l')

# R
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[6,], type='l', xlim=rev(range(brNesovin$t)))
lines(br$t + br$depth, br$asr[6,], type='l')


513-Ntip(phy)
495-Ntip(phy)
st[,242]
st[,224]

# --- Plot setup ---
# states <- c("A", "E", "M", "U", "S", "R",
#             "A.E", "A.M", "A.U", "A.S", "A.R",
#             "E.M", "E.U", "S.R")
# plot(NULL,
#      xlim = rev(range(br$t + br$depth)),
#      ylim = c(0, 1),
#      xlab = "Time before present (Ma)",
#      ylab = "Probability",
#      main = "Ancestral State Probabilities through Time",
#      las = 1, cex.lab = 1.2, cex.axis = 1)
#
# # --- Choose focal states to plot ---
# focal_states <- c("A", "S", "R", "A.S", "S.R")
# focal_indices <- match(focal_states, states)
#
# # --- Add lines for selected states ---
# for (i in focal_indices) {
#   lines(br$t + br$depth, br$asr[i, ],
#         col = colsDB[states[i]],
#         lwd = 2)
# }
#
# # --- Add legend ---
# legend("topright", legend = focal_states,
#        col = colsDB[focal_states],
#        lty = 1, lwd = 2, bty = "n", cex = 0.9)

#--------------
focal_states <- c("A", "A.S", "S", "S.R", "R")
focal_indices <- match(focal_states, states)

# Set up the plotting area: 5 rows, 1 column
par(mfcol = c(5,1), mar = c(1.5, 6, 1, 5), oma = c(4, 0, 2, 0))  # smaller vertical

#par(mfcol = c(5,1), mar = c(3.5, 4, 1, 1), oma = c(0, 0, 2, 0))  # adjust margins

for (idx in seq_along(focal_indices)) {
  i <- focal_indices[idx]
  x <- br$t + br$depth
  y <- br$asr[i, ]

  if (idx != length(focal_indices)) {
    plot(x, y, type = "n",
         xaxt = "n", yaxt = "n",
         xlab = NA,
         ylab = paste("", states[i]),
         xlim = rev(range(x)),
         ylim = range(br$asr[i, ]))
    #axis(side = 1, labels = FALSE)  # x-ticks only
  } else {
    plot(x, y, type = "n",
         xlab = "Time before present (Ma)",
         yaxt = "n",
         ylab = paste("", states[i]),
         xlim = rev(range(x)),
         ylim = range(br$asr[i, ]))
  }

  # Y-axis: show only 0 and max
  axis(side = 2, at = c(0, max(y)/2, max(y)), labels = c("0", NA, round(max(y), 2)) )

  polygon(c(x, rev(x)), c(rep(0, length(y)), rev(y)),
          col = adjustcolor(colsDB[states[i]], alpha.f = 0.3), border = NA)
  lines(x, y, col = colsDB[states[i]], lwd = 2)
}



#----------- Nesovinsonia

#---------------------------------
#----------- Nesovinsonia
extract.clade.simm<-function(tree,node){
  x<-getDescendants(tree,node)
  x<-x[x<=length(tree$tip.label)]
  drop.tip(tree,tree$tip.label[-x])
}
# ASR plot
x<-getDescendants(phy, 489)
#Descendants(phy, 495, type = c( "children"))
x<-x[x >length(phy$tip.label)]
length(x)
x <- c(489, x)
x <- x[order(x)]
#Descendants(x, node, type = c("tips", "children", "all"))
phySis <- extract.clade.simm(phy, 489)
stSis <- st[,x-Ntip(phy)]

geo_regionsSis <- geo_regions[match(phySis$tip.label, names(geo_regions))]

par(mfcol = c(5,1), mar = c(3, 2, 3, 3), oma = c(3, 3, 2, 1))
pdf("runs/asr/asr-Nesovinsonia.pdf", width = 10, height = 3)  # Change size as needed
plotTree(phySis, fsize=1.5, ftype="i", lwd=5, offset=.8)
nodelabels(node=1:phySis$Nnode+Ntip(phySis), pie=t(stSis), piecol=colsDB, cex=1)
#nodelabels(cex = .3, height = .05, width = .05, frame='none', col ='grey')
tiplabels(pie=to.matrix(geo_regionsSis,sort(unique(geo_regionsSis))), piecol=colsDB[-c(2,3,4,5)], cex=0.6, adj = 0.5)
#axisPhylo()
# --- Close device ---
dev.off()

phy$tip.label
# 489  Neso + Tanzanolus

brNesovin <- asr.marginal.classe_branch.multiple(lik.td, mle.td$par.full, node.id=67, eps = 0.1, Nbins=100, include.ends=TRUE)

# states: [A, E, M, U, S, R,   A.E, A.M, A.U, A.S, A.R,   E.M, E.U,    S.R]

# A
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[1,], type='l', xlim=rev(range(brNesovin$t)))
# A.S
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[10,], type='l', xlim=rev(range(brNesovin$t)))
# S
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[5,], type='l', xlim=rev(range(brNesovin$t)))
# S.R
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[14,], type='l', xlim=rev(range(brNesovin$t)))
# R
plot(brNesovin$t + brNesovin$depth, brNesovin$asr[6,], type='l', xlim=rev(range(brNesovin$t)))

#--------------
focal_states <- c("A", "A.S", "S", "S.R", "R")
focal_indices <- match(focal_states, states)

# Set up the plotting area: 5 rows, 1 column
par(mfcol = c(5,1), mar = c(1.5, 6, 1, 5), oma = c(4, 0, 2, 0))  # smaller vertical

#par(mfcol = c(5,1), mar = c(3.5, 4, 1, 1), oma = c(0, 0, 2, 0))  # adjust margins

for (idx in seq_along(focal_indices)) {
  i <- focal_indices[idx]
  x <- brNesovin$t + brNesovin$depth
  y <- brNesovin$asr[i, ]

  if (idx != length(focal_indices)) {
    plot(x, y, type = "n",
         xaxt = "n", yaxt = "n",
         xlab = NA,
         ylab = paste("", states[i]),
         xlim = rev(range(x)),
         ylim = range(brNesovin$asr[i, ]))
    #axis(side = 1, labels = FALSE)  # x-ticks only
  } else {
    plot(x, y, type = "n",
         xlab = "Time before present (Ma)",
         yaxt = "n",
         ylab = paste("", states[i]),
         xlim = rev(range(x)),
         ylim = range(brNesovin$asr[i, ]))
  }

  # Y-axis: show only 0 and max
  axis(side = 2, at = c(0, max(y)/2, max(y)), labels = c("0", NA, round(max(y), 2)) )

  polygon(c(x, rev(x)), c(rep(0, length(y)), rev(y)),
          col = adjustcolor(colsDB[states[i]], alpha.f = 0.3), border = NA)
  lines(x, y, col = colsDB[states[i]], lwd = 2)
}



#---- Expected N of species
lam <- 0.232526
mu <- 0.473558
t <- 15

1*exp(lam-mu)*5

1*exp(lam-mu)*40





#---- Mada

# Epacto
brEpa<- asr.marginal.classe_branch.multiple(lik.td, mle.td$par.full, node.id=520, eps = 0.1, Nbins=100, include.ends=TRUE)
brNano<- asr.marginal.classe_branch.multiple(lik.td, mle.td$par.full, node.id=443, eps = 0.1, Nbins=100, include.ends=TRUE)
brEpi<- asr.marginal.classe_branch.multiple(lik.td, mle.td$par.full, node.id=394, eps = 0.1, Nbins=100, include.ends=TRUE)
brHel<- asr.marginal.classe_branch.multiple(lik.td, mle.td$par.full, node.id=291, eps = 0.1, Nbins=100, include.ends=TRUE)


# states: [A, E, M, U, S, R,   A.E, A.M, A.U, A.S, A.R,   E.M, E.U,    S.R]

par(mfcol = c(5,1), mar = c(1.5, 6, 1, 5), oma = c(4, 0, 2, 0))
plot(brEpa$t + brEpa$depth, brEpa$asr[3,], type='l', xlim=rev(range(brEpa$t+ brEpa$depth)))
plot(brNano$t + brNano$depth, brNano$asr[3,], type='l', xlim=rev(range(brNano$t+ brNano$depth)))
plot(brEpi$t + brEpi$depth, brEpi$asr[3,], type='l', xlim=rev(range(brEpi$t+ brEpi$depth)))
plot(brHel$t + brHel$depth, brHel$asr[3,], type='l', xlim=rev(range(brHel$t+ brHel$depth)))

par(mfcol = c(5,1), mar = c(1.5, 6, 1, 5), oma = c(4, 0, 2, 0))
plot(brEpa$t + brEpa$depth, brEpa$asr[3,], type='l', xlim=c(20, 5))
plot(brNano$t + brNano$depth, brNano$asr[3,], type='l', xlim=c(20, 5))
plot(brEpi$t + brEpi$depth, brEpi$asr[3,], type='l', xlim=c(20, 5))
plot(brHel$t + brHel$depth, brHel$asr[3,], type='l', xlim=c(20, 5))

# Set up plotting area
par(mfcol = c(5, 1), mar = c(1.5, 6, 1, 5), oma = c(4, 0, 2, 0))

# Define plot data list
plot_data <- list(
  list(data = brEpa, label = "Epa"),
  list(data = brNano, label = "Nano"),
  list(data = brEpi, label = "Epi"),
  list(data = brHel, label = "Hel")
)

# Use uniform color
uniform_col <- "#74ADD1"

# Loop through datasets
for (idx in seq_along(plot_data)) {
  br <- plot_data[[idx]]$data
  label <- plot_data[[idx]]$label
  x <- br$t + br$depth
  y <- br$asr[3, ]

  # Set up empty plot
  plot(x, y, type = "n",
       xlim = c(20, 5),  # reversed axis
       ylim = range(y),
       xaxt = if (idx != length(plot_data)) "n" else "s",
       yaxt = "n",
       xlab = if (idx == length(plot_data)) "Time before present (Ma)" else NA,
       ylab = NA)
  axis(side = 1, labels = FALSE)

  # Add filled polygon
  polygon(c(x, rev(x)), c(rep(0, length(y)), rev(y)),
          col = adjustcolor(uniform_col, alpha.f = 0.3), border = NA)

  # Add line on top
  lines(x, y, col = uniform_col, lwd = 2)

  # Add custom y-axis: only 0 and max
  axis(side = 2,
       at = c(min(y),   max(y)),
       labels = c(round(min(y), 2),   round(max(y), 2)),
       las = 1)  # vertical tick labels
}
