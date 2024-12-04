
comb2pars <- function(combReg, args){
  lamT <- args$lam.tensor
  parameters <- c()
  for (i in seq_along(combReg)){
    #print(combReg[i])
    par_values <- combReg[i][[1]]
    reg_name=names(combReg[i])
    split_string <- strsplit(reg_name, ".", fixed = TRUE)[[1]]
    y = split_string[1]
    x = split_string[2]
    par_names <- c(
      lamT[[reg_name]][y,x],
      lamT[[reg_name]][y,reg_name],
      lamT[[reg_name]][x,reg_name]
    )
    names(par_values) <- par_names
    parameters <- c(parameters, par_values)
  }
  return(parameters)
}


#unitReg <- yaml$lambda_unit_regions
unit2pars <- function(unitReg, args){
  lamT <- args$lam.tensor
  parameters <- c()
  for (i in seq_along(unitReg)){
    #print(combReg[i])
    par_values <- unitReg[i][[1]]
    reg_name=names(unitReg[i])
    par_names <- lamT[[reg_name]][reg_name,reg_name]
    names(par_values) <- par_names
    parameters <- c(parameters, par_values)
  }
  return(parameters)
}

mu2pars <- function(unitReg, args){
  mu <- args$mu
  parameters <- c()
  i=1
  for (i in seq_along(unitReg)){
    #print(combReg[i])
    par_values <- unitReg[i][[1]]
    reg_name=names(unitReg[i])
    par_names <- mu[reg_name]
    names(par_values) <- par_names
    parameters <- c(parameters, par_values)
  }
  return(parameters)
}


# Function to remove "0" elements from a vector
remove_zeros <- function(vec) {
  vec[vec != "0"]
}


dQ2pars <- function(dQ, args){
  Q <- args$Q
  indices <- !(colnames(Q) %in% names(dQ))
  to_reg_names <- colnames(Q)[indices]

  parameters <- c()
  i=1
  for (i in seq_along(dQ)){
    #print(dQ[i])
    par_values <- dQ[i][[1]]
    reg_name=names(dQ[i])
    #split_string <- strsplit(reg_name, ".", fixed = TRUE)[[1]]
    par_names <- Q[reg_name,to_reg_names]
    names(par_values) <- par_names
    par_values <- remove_zeros(par_values)
    parameters <- c(parameters, par_values)
  }
  return(parameters)
}

#yaml_file="my-test/data/2-regions.yml"
#' Title
#'
#' @param yaml_file
#' @param args
#'
#' @return
#' @export
#'
#' @examples
yaml2vec <- function(yaml_file, args){

  to_char <- function(x){as.character(x)}
  yaml <- yaml.load_file(yaml_file, handlers=list("int"= to_char))
  #yaml
  unitReg <- yaml$lambda_unit_regions
  sA <- unit2pars(unitReg, args)

  combReg <- yaml$lambda_comb_regions
  sAB <- comb2pars(combReg, args)

  dQ <- yaml$Q_dQ
  dX <- dQ2pars(dQ, args)

  unitReg <- yaml$mu
  mu <- mu2pars(unitReg, args)

  xQ <- yaml$Q_xQ
  xX <- dQ2pars(xQ, args)

  par.vec <- c(sA, sAB, dX, mu, xX)
  grouped_pars <- split(names(par.vec ), par.vec )

  out <- list(Nstates=length(yaml$states), states=yaml$states, pars=grouped_pars)

  #return(c(sA, sAB, dX, mu, xX))
  return(out)
}

# yaml_file="my-test/data/2-regions.yml"
# read_yaml_pars(yaml_file)
read_yaml_pars <-  function(yaml_file){

  to_char <- function(x){as.character(x)}
  yaml <- yaml.load_file(yaml_file, handlers=list("int"= to_char))

  regions <- yaml$states
  Nstates <- length(yaml$states)
  pars <- diversitree:::default.argnames.classe(Nstates)
  args <- pars_to_arrays(pars, Nstates, regions)
  #args
  yaml2vec(yaml_file, args)
}

# yaml_file="my-test/data/2-regions-td.yml"
# read_yaml_pars(yaml_file)
read_yaml_pars_td <-  function(yaml_file){

  # read yaml
  to_char <- function(x){as.character(x)}
  yaml <- yaml.load_file(yaml_file, handlers=list("int"= to_char))

  # get metadata
  regions <- yaml$states
  Nstates <- length(yaml$states)
  n.epoch <- as.numeric(yaml$n.epoch)
  epoch.times <- as.numeric(yaml$epoch.times)

  #pars <- diversitree:::default.argnames.classe(Nstates)
  #args <- pars_to_arrays(pars, Nstates, regions)
  #args
  # argnames_classe = diversitree:::default.argnames.classe(Nstates) #length(argnames_classe)
  # argnames.td <- c(sprintf("t.%d", seq_len(n.epoch - 1)), diversitree:::argnames_twopart(argnames_classe, n.epoch))

  # genertae args
  args <- vector("list", length = n.epoch)
  for (epoch in seq_along(1:n.epoch)){
    pars <- diversitree:::default.argnames.classe(Nstates)
    pars <- paste0(pars, '.', epoch)
    args.i <- pars_to_arrays(pars, Nstates, regions)
    args[[epoch]] <- args.i
  }

  #yaml parameters
  # epoch=1
  grouped_pars <- vector("list", length = n.epoch)
  for (epoch in seq_along(1:n.epoch)){
    yaml_epoch_i <- yaml[[as.character(epoch)]]

    unitReg <- yaml_epoch_i$lambda_unit_regions
    sA <- unit2pars(unitReg, args[[epoch]])

    combReg <- yaml_epoch_i$lambda_comb_regions
    sAB <- comb2pars(combReg, args[[epoch]])

    dQ <- yaml_epoch_i$Q_dQ
    dX <- dQ2pars(dQ, args[[epoch]])

    unitReg <- yaml_epoch_i$mu
    mu <- mu2pars(unitReg, args[[epoch]])

    xQ <- yaml_epoch_i$Q_xQ
    xX <- dQ2pars(xQ, args[[epoch]])

    par.vec <- c(sA, sAB, dX, mu, xX)
    grouped_pars_i <- split(names(par.vec ), par.vec )
    grouped_pars[[epoch]] <- grouped_pars_i
  }

  # Combine elements with the same names across sublists
  combined_list <- lapply(unique(unlist(lapply(grouped_pars, names))), function(name) {
    unique(unlist(lapply(grouped_pars, function(x) x[[name]])))
  })

  # Assign names to the combined list
  names(combined_list) <- unique(unlist(lapply(grouped_pars, names)))

  out <- list(Nstates=length(yaml$states), states=yaml$states, n.epoch=n.epoch, epoch.times=epoch.times,  pars=combined_list)


  return(out)
}


# yaml_file="my-test/data/2-regions.yml"
# par.categories <- read_yaml_pars(yaml_file)
# make_constraints_sse(par.categories)
make_constraints_sse <- function(par.categories){
  #regions <- yaml$states
  Nstates <- par.categories$Nstates
  pars <- diversitree:::default.argnames.classe(Nstates)
  pars2formula(par.categories$pars, pars)
}


# par.categories.td <- read_yaml_pars_td("my-test/data/2-regions-td.yml")
# make_constraints_sse_td(par.categories)
make_constraints_sse_td <- function(par.categories.td){
  #regions <- yaml$states
  Nstates <- par.categories.td$Nstates
  pars <- diversitree:::default.argnames.classe(Nstates)
  pars2formula(par.categories$pars, pars)
  #------
  argnames_classe <-  diversitree:::default.argnames.classe(Nstates) #length(argnames_classe)
  pars <- diversitree:::argnames_twopart(argnames_classe, n.epoch)
  # argnames.td <- c(sprintf("t.%d", seq_len(n.epoch - 1)), diversitree:::argnames_twopart(argnames_classe, n.epoch))

  # Formulas for the times of epochs
  times <- par.categories.td$epoch.times
  parameters <- sprintf("t.%d", seq_len(n.epoch - 1))
    # Create the formulas
  time_formulas <- mapply(function(param, time) as.formula(paste(param, "~", time)),  parameters, times, SIMPLIFY = FALSE)

  # other formulas
  existing_formulas <- pars2formula(par.categories.td$pars, pars)

  # Append the new formulas to the existing ones
  all_formulas <- c(existing_formulas, new_formulas)

  return(all_formulas)

}

# reg <- c("A", "B", "A.B")
# pars <- diversitree:::default.argnames.classe(3)
# args <- pars_to_arrays(pars, 3, reg)
# grouped_pars <- yaml2vec("my-test/data/2-regions.yml", args)
# view_pars_as_array(pars, grouped_pars, k=3, regions=reg)

view_pars_as_array <- function(par.categories){

  k <- par.categories$Nstates
  regions <- par.categories$states
  pars <- diversitree:::default.argnames.classe(k)
  grouped_pars <- par.categories$pars

  pars.constr <- rep(0, length(pars))
  names(pars.constr) <- pars

  for (i in 1:length(grouped_pars)){
    par.group <- grouped_pars[[i]]
    par.group.name <- names(grouped_pars[i])
    for (j in 1:length(par.group)){
      #print(par.group[j])
      pars.constr[[par.group[j]]] <-  par.group.name
    }
  }
  pars_to_arrays(pars.constr, k, regions)

}


#------------------------
get_single_pars <- function(list_of_vectors){
  # Select vectors with length 1
  single_length_vectors <- lapply(list_of_vectors, function(x) {
    if(length(x) == 1) {
      return(x)
    }
  })
  # Remove NULL elements (vectors with length != 1)
  single_length_vectors <- single_length_vectors[!sapply(single_length_vectors, is.null)]
  # Output the list of vectors with length 1
  unlist(single_length_vectors)
}

get_multi_pars <- function(list_of_vectors){
  # Select vectors with length 1
  single_length_vectors <- lapply(list_of_vectors, function(x) {
    if(length(x) > 1) {
      return(x)
    }
  })
  # Remove NULL elements (vectors with length != 1)
  single_length_vectors <- single_length_vectors[!sapply(single_length_vectors, is.null)]
  # Output the list of vectors with length 1
  single_length_vectors
}


create_formulas <- function(vec_list) {
  all_formulas <- list()
  for (vec in vec_list) {
    if (length(vec) > 1) {
      for (i in 2:length(vec)) {
        formula <- paste(vec[i], "~", vec[1])
        all_formulas <- c(all_formulas, list(formula))
      }
    }
  }
  return(all_formulas)
}

# grouped_pars <- yaml2vec("my-test/data/2-regions.yml", args)
#' Title
#'
#' @param grouped_pars
#' @param pars
#'
#' @return
#' @export
#'
#' @examples
pars2formula <- function(grouped_pars, pars){
  #grouped_pars <- split(names(par.vec ), par.vec )
  #single_pars <- get_single_pars(grouped_pars)
  multi_pars <- get_multi_pars(grouped_pars)
  formula.multi <- create_formulas(multi_pars)

  non_zero <-unlist(grouped_pars)
  #non_zero <- terms_from_formula(f.list)
  zero_terms <- pars[!pars %in% non_zero]
  formula.zero <-lapply(zero_terms, function(arg) reformulate("0", response = arg))
  formula.out <- c(formula.zero, formula.multi)

  return(formula.out)
}
