



par.categories.td <- read_yaml_pars_td("my-test/yml/geosse-4reg-BC.yml")
par.symbolic.td <- pars_yaml_to_arrays_td(par.categories.td)
print(par.categories.td)
print(par.symbolic.td)

create_parameter_html(par.symbolic.td[[1]])


reg <- c("A", "B", "C", "A.B", "A.C", "B.C")
pars <- diversitree:::default.argnames.classe(6)
args <- pars_to_arrays(pars, 6, reg)
print(args)


source('R/html.R')
#create_parameter_html(par.symbolic.td[[1]])
create_parameter_html(par.symbolic.td[[1]], cell_names = args)


#--- 2 epochs
source('R/html.R')
#create_parameter_html(par.symbolic.td[[1]])
par.categories.td <- read_yaml_pars_td("my-test/yml/geosse-4reg-BC.yml")
par.symbolic.td <- pars_yaml_to_arrays_td(par.categories.td)
create_parameter_html(par.symbolic.td, cell_names = list(args, args), output_file = "model_parameters1.html")

source('R/html.R')
par.categories.td <- read_yaml_pars_td("my-test/yml/geosse-4reg-BC.yml")
create_parameter_html(par.categories.td, "mode1", output_file = "model_parameters2.html")


Nstates <- par.categories.td$Nstates
n.epoch <- par.categories.td$n.epoch
epoch.times <- par.categories.td$epoch.times
regions <-par.categories.td$states
length(par.categories.td$pars)

# genertae parameter names
cell_names <- vector("list", length = n.epoch)
for (epoch in seq_along(1:n.epoch)){
  pars <- diversitree:::default.argnames.classe(Nstates)
  pars <- paste0(pars, '.', epoch)
  cell_names.i <- pars_to_arrays(pars, Nstates, regions)
  cell_names[[epoch]] <- cell_names.i
}


#-- BC
par.cat4 <- read_yaml_pars_td("my-test/yml/geosse-4reg-BC.yml")
print(par.cat4)


reg <- c("A", "B", "C", "A.B", "A.C", "B.C")
pars <- diversitree:::default.argnames.classe(6)
args <- pars_to_arrays(pars, 6, reg)
print(args)

create_parameter_html(args)

#-- CB

par.cat4.CB <- read_yaml_pars_td("yml/geosse-4reg-CB.yml")
print(par.cat4.CB)


reg <- c("A", "B", "C", "A.B", "A.C", "C.B")
pars <- diversitree:::default.argnames.classe(6)
args.CB <- pars_to_arrays(pars, 6, reg)
print(args.CB)

