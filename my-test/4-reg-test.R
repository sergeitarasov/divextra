
#-- BC

par.cat4 <- read_yaml_pars_td("yml/geosse-4reg-BC.yml")
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

