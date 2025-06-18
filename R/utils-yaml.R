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
yaml2vec <- function(yaml_file, args){

  to_char <- function(x){as.character(x)}
  yaml <- yaml::yaml.load_file(yaml_file, handlers=list("int"= to_char))
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


#' Read parameters of time-homogeneous GeoSSE model from yaml file
#'
#' @param yaml_file yaml file
#'
#' @return named list of parameters and configurations
#' @export
#'
#' @examples
#'
#' file_path <- system.file("extdata", "geosse3_ti_maxpar.yml", package = "divextra")
#' par.categories <- read_yaml_pars(file_path)
#' print(par.categories)
read_yaml_pars <-  function(yaml_file){

  to_char <- function(x){as.character(x)}
  yaml <- yaml::yaml.load_file(yaml_file, handlers=list("int"= to_char))

  regions <- yaml$states
  Nstates <- length(yaml$states)
  pars <- diversitree:::default.argnames.classe(Nstates)
  args <- pars_to_arrays(pars, Nstates, regions)
  #args
  yaml2vec(yaml_file, args)
}


check_zero_in_list <- function(lst) {
  # Recursive function to check all elements
  check_element <- function(x) {
    if (is.list(x)) {
      # If element is a list, recursively check its elements
      lapply(x, check_element)
    } else {
      # Convert element to character and check for exact "0"
      if (any(x == "0")) {
        print(lst)
        stop("Error: List contains standalone '0' element.\n
        Revise the order of dual states in yaml.\n
        Correct:   A, B, C, A.B, A.C,   B.C\n
        Incorrect: A, B, C, A.B, A.C,   C.B\n")
      }
    }
  }

  # Check all elements in the main list
  invisible(lapply(lst, check_element))

  # Return TRUE if no standalone zeros found
  return(TRUE)
}


# check_zero_in_list <- function(lst) {
#   # Recursive function to check all elements
#   check_element <- function(x) {
#     if (is.list(x)) {
#       # If element is a list, recursively check its elements
#       lapply(x, check_element)
#     } else {
#       # Convert element to character and check for "0"
#       if (any(grepl("0", as.character(x)))) {
#         print(lst)
#         stop("Error: List contains element with '0' character.\n
#         Revise the order of dual states in yaml.\n
#         Correct:   A, B, C, A.B, A.C,   B.C\n
#         Incorrect: A, B, C, A.B, A.C,   C.B\n")
#       }
#     }
#   }

#   # Check all elements in the main list
#   invisible(lapply(lst, check_element))

#   # Return TRUE if no zeros found
#   return(TRUE)
# }

#' Read parameters of episodic GeoSSE model from yaml file
#'
#' @param yaml_file yaml file
#'
#' @return named list of parameters and configurations
#' @export
#'
#' @examples
#'
#' file_path <- system.file("extdata", "geosse3_td_maxpar.yml", package = "divextra")
#' par.categories.td <- read_yaml_pars_td(file_path)
#' print(par.categories.td)
read_yaml_pars_td <-  function(yaml_file){

  # read yaml
  to_char <- function(x){as.character(x)}
  yaml <- yaml::yaml.load_file(yaml_file, handlers=list("int"= to_char))

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

  # remove zero valued pars: they should be automatically set to 0's in the ln constrain function
  combined_list$`0` <- NULL
  # check if any pars are zero values (they should not!
  # but it may happen if the state oreder is incorrect in yaml)
  check_zero_in_list(combined_list)

  out <- list(Nstates=length(yaml$states), states=yaml$states, n.epoch=n.epoch, epoch.times=epoch.times,  pars=combined_list)


  return(out)
}


#' Create parameter constraints for time-homogeneous ClaSSE
#'
#' @param par.categoriesa a named list of parameters, for example, from read_yaml_pars_td()
#'
#' @return list of formulas
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "geosse3_ti_maxpar.yml", package = "divextra")
#' par.categories <- read_yaml_pars(file_path)
#' print(par.categories)
#' pars <- make_constraints_sse(par.categories)
#' print(pars)
make_constraints_sse <- function(par.categories){
  #regions <- yaml$states
  Nstates <- par.categories$Nstates
  pars <- diversitree:::default.argnames.classe(Nstates)
  pars2formula(par.categories$pars, pars)
}





#' Create parameter constraints for episodic ClaSSE
#'
#' @param par.categories.td a named list of parameters, for example, from read_yaml_pars_td()
#'
#' @return list of formulas
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "geosse3_td_maxpar.yml", package = "divextra")
#' par.categories.td <- read_yaml_pars_td(file_path)
#' print(par.categories.td)
#' formula.td <- make_constraints_sse_td(par.categories.td)
#' print(formula.td)
make_constraints_sse_td <- function(par.categories.td){
  #regions <- yaml$states
  n.epoch = par.categories.td$n.epoch
  Nstates <- par.categories.td$Nstates
  #pars <- diversitree:::default.argnames.classe(Nstates)
  #pars2formula(par.categories$pars, pars)

  #
  argnames_classe <-  diversitree:::default.argnames.classe(Nstates) #length(argnames_classe)
  pars <- diversitree:::argnames_twopart(argnames_classe, n.epoch)
  # argnames.td <- c(sprintf("t.%d", seq_len(n.epoch - 1)), diversitree:::argnames_twopart(argnames_classe, n.epoch))

  # Formulas for the times of epochs
  times <- par.categories.td$epoch.times
  parameters <- sprintf("t.%d", seq_len(n.epoch - 1))
    # Create the formulas
  time_formula <- mapply(function(param, time) as.formula(paste(param, "~", time)),  parameters, times, SIMPLIFY = FALSE)

  # other formulas
  existing_formulas <- pars2formula(par.categories.td$pars, pars)

  # Append the new formulas to the existing ones
  all_formulas <- c(existing_formulas, time_formula)

  return(all_formulas)

}



#' Display symbolic parameters of time-homogeneous ClaSSE as arrays
#'
#' @param par.categories named list of parameters and configurations (e.g., from yaml file)
#'
#' @return a list of arrays
#' @export
#'
#' @examples
#'
#' file_path <- system.file("extdata", "geosse3_ti_maxpar.yml", package = "divextra")
#' yaml_pars <- read_yaml_pars(file_path)
#' pars_yaml_to_arrays(yaml_pars)
pars_yaml_to_arrays <- function(par.categories){

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



# Function to filter elements in the list
filter_elements_by_suffix <- function(named_list, suffix) {
  lapply(named_list, function(vec) vec[grepl(paste0("\\", suffix, "$"), vec)])
}

# Function to remove tags from elements in the list
remove_suffix <- function(named_list, suffix) {
  lapply(named_list, function(vec) gsub(paste0("\\", suffix, "$"), "", vec))
}

# Function to remove list elements that are character(0)
remove_empty_list_elements <- function(lst) {
  lst[sapply(lst, function(x) length(x) > 0)]
}


#' Display symbolic parameters of episodic ClaSSE as arrays
#'
#' @param par.categories.td named list of parameters and configurations (e.g., from yaml file)
#'
#' @return a list of arrays
#' @export
#'
#' @examples
#'
#' file_path <- system.file("extdata", "geosse3_td_test.yml", package = "divextra")
#' par.categories.td <- read_yaml_pars_td(file_path)
#' par.symbolic.td <- pars_yaml_to_arrays_td(par.categories.td)
#' print(par.categories.td)
#' print(par.symbolic.td)
pars_yaml_to_arrays_td <- function(par.categories.td){
  new.pars <- par.categories.td
  n.epochs <- par.categories.td$n.epoch
  out <- vector('list', n.epochs)
  # i=1
  for (i in 1:n.epochs){
    suffix <- paste0('.', i)
    # Remove elements not ending with e.g. '.1'
    filtered_list <- filter_elements_by_suffix(par.categories.td$pars, suffix)
    # Remove e.g., '.1' from elements
    xxx <- remove_suffix(filtered_list, suffix)
    # remove empty sublists
    new.pars$pars <-remove_empty_list_elements(xxx)
    out[[i]] <- pars_yaml_to_arrays(new.pars)
  }
  return(out)
}



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





are_lists_equivalent <- function(list1, list2) {
  # Helper function to standardize a vector
  standardize_vector <- function(vec) {
    if (length(vec) > 0) {
      return(paste(sort(vec), collapse = "|"))
    }
    return("")
  }

  # Convert each vector to a sorted string representation
  strings1 <- sapply(list1, standardize_vector)
  strings2 <- sapply(list2, standardize_vector)

  # Sort the string representations
  strings1 <- sort(strings1)
  strings2 <- sort(strings2)

  # Compare the sorted strings
  return(identical(strings1, strings2))
}

# Example usage:
# l1 <- list(c("a1", "bc", "cwe"), c("a"))
# l2 <- list(c("a"), c("bc", "a1", "cwe"))
# are_lists_equivalent(l1, l2)  # Returns TRUE
#
# l3 <- list(c("a1", "bc", "cwe"), c("a", "d"))
# are_lists_equivalent(l1, l3)  # Returns FALSE

# Function to remove elements containing "t." from a list
remove_t_elements <- function(lst) {
  # Find which elements contain "t."
  has_t <- sapply(lst, function(x) any(grepl("t\\.", x)))

  # Return list without those elements
  lst[!has_t]
}

# Example usage:
# my_list <- list(
#   "0.655" = c("mu01.1", "mu02.1"),
#   "0.703" = c("lambda010101.1", "lambda020202.1"),
#   "10" = c("t.1")
# )
# cleaned_list <- remove_t_elements(my_list)




#' Check consistency between diversitree and yaml parameter groups
#'
#' @description
#' Compares parameter groupings from diversitree MLE results with parameter groups defined in yaml configuration.
#' Removes zero-valued parameters and timing parameters before comparison.
#'
#' @param par.diversitree List of parameter groups from diversitree MLE results
#' @param par.yaml List of parameter groups from yaml configuration
#'
#' @return Logical. TRUE if parameter groupings are equivalent, FALSE otherwise
#' @export
#'
#' @details
#' The function performs the following steps:
#' 1. Removes zero-valued parameters ('0' group)
#' 2. Removes timing parameters (containing 't.')
#' 3. Removes list names for comparison
#' 4. Checks if remaining parameter groups are equivalent regardless of order
#'
#' @examples
#' \dontrun{
#' par.categories.td <- read_yaml_pars_td("yml/geosse.yml")
#' mle.td <- find.mle(lik.const.td, starting.point, condition.surv=TRUE)
#' mle.pars.groups <- pars2groups(mle.td$par.full)
#' check_parametrization_consistency(mle.pars.groups, par.categories.td$pars)
#' }
check_parametrization_consistency <- function(par.diversitree, par.yaml){
  par.diversitree$'0' <- NULL
  par.diversitree <- remove_t_elements(par.diversitree)
  names(par.diversitree)<-NULL
  names(par.yaml) <- NULL

  are_lists_equivalent(par.diversitree, par.yaml)
}
