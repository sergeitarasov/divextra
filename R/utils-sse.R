
fill_off_diagonal <- function(off_diagonal_elements, N) {
  # Determine N based on length of off_diagonal_elements
  #N <- (1 + sqrt(1 + 8 * length(off_diagonal_elements))) / 2
  #N <- as.integer(N)

  # Create an N x N zero matrix
  matrix <- matrix(0, nrow=N, ncol=N)

  counter <- 1

  for (i in 1:N) {
    for (j in 1:N) {
      if (i != j) {
        matrix[i, j] <- off_diagonal_elements[counter]
        counter <- counter + 1
      }
    }
  }

  return(matrix)
}
# Test
# QQ <- generateRateMatrix(3)
# QQ
# extr <- extract_off_diagonal(QQ)
# extr
# fill_off_diagonal(extr, 3)


generateUpperTriangularMatrix_classe <- function(element, N) {
  # Create an empty matrix
  matrix_result <- matrix(numeric(0), nrow = N, ncol = N)

  # Populate the matrix with indexed elements or 0 based on position
  # i=1
  # j=2
  for (i in 1:N) {
    for (j in 1:N) {
      if (i <= j) {
        if (class(element)=="character")
          matrix_result[i, j] <- paste0(element, i, j)
        else
          matrix_result[i, j] <- (element)
      } else {
        matrix_result[i, j] <- 0
      }
    }
  }

  return(matrix_result)
}
# # Test
# generateUpperTriangularMatrix_classe('lam0', 3)
# generateUpperTriangularMatrix_classe(1, 3)



vector_to_upper_triangular_matrix <- function(N) {
  # Compute the side length of the matrix
  n <- (-1 + sqrt(1 + 8 * length(N))) / 2
  n <- as.integer(n)

  # Initialize a zero matrix of size n x n
  matrix <- matrix(0, n, n)

  counter <- 1
  for(i in 1:n) {
    for(j in i:n) {
      matrix[i, j] <- N[counter]
      counter <- counter + 1
    }
  }

  return(matrix)
}
# Test
# vector_to_upper_triangular_matrix(c(1,2,3,4,5,6))



split_vector <- function(vec, K) {
  # Create a sequence repeating from 1 to (N/K) for each element of vec
  f <- rep(1:(length(vec)/K), each=K)

  # Split the vector based on the sequence
  list_of_vectors <- split(vec, f)

  return(list_of_vectors)
}
# # Test
# vec <- 1:20  # A vector with 20 elements
# K <- 5       # We want to split it into vectors of size 5
# result <- split_vector(vec, K)
# result





#' Display ClaSSE parameters as arrays
#'
#' @description
#' Takes a named vector of ClaSSE parameters and displays them as arrays
#'
#'
#' @param pars ClaSSE parameters
#' @param Nstates number of character states
#' @param nm names for character states
#'
#' @return a list of arrays
#' @export
#'
#' @examples
#'
#' reg <- c("A", "B", "A.B")
#' pars <- diversitree:::default.argnames.classe(3)
#' args <- pars_to_arrays(pars, 3, reg)
#' print(args)
pars_to_arrays <- function(pars, Nstates, nm=NULL){
  # number of lambdas
  #0.5*(Nstates*Nstates-Nstates)+Nstates
  Nlambdas_one_state <- (0.5*(Nstates*Nstates+Nstates))
  Nlambdas <- Nlambdas_one_state*Nstates
  #Nqs <- Nstates*Nstates - Nstates

  lambdas <- pars[1:Nlambdas]
  mu <- pars[(Nlambdas+1):(Nlambdas+Nstates)]
  qs <- pars[-1:-(Nlambdas+Nstates)]

  lam.tensor <- split_vector(lambdas, Nlambdas_one_state)
  lam.tensor <- lapply(lam.tensor, function(x) vector_to_upper_triangular_matrix(x))
  Q <- fill_off_diagonal(qs, Nstates)

  #if (is.numeric(Q[1,1])) diag(Q) <- -rowSums(Q)

  # add names
  # nm <- c('a','b','c')
  # Function to set column and row names
  set_names <- function(mat, names) {
    colnames(mat) <- names
    rownames(mat) <- names
    return(mat)
  }
  lam.tensor <- lapply(lam.tensor, set_names, names = nm)
  names(lam.tensor) <- nm
  names(mu) <- nm
  rownames(Q) <- colnames(Q) <- nm

  return(list(lam.tensor=lam.tensor, mu=mu, Q=Q))
}



#' Reorder character states in arrays
#'
#' @param pars.array a list of arrays
#' @param v new order of character states
#'
#' @return a list of arrays
#' @export
#'
#' @examples
#'
#' reg <- c("A", "B", "A.B")
#' pars <- diversitree:::default.argnames.classe(3)
#' args <- pars_to_arrays(pars, 3, reg)
#' print(args)
#' new.order <- reoder_pars(args, c(3,1,2))
#' print(new.order)
reoder_pars <- function(pars.array, v){
  lam.tensor <- pars.array$lam.tensor
  mu <- pars.array$mu
  Q <- pars.array$Q

  lam.tensor <- lapply(lam.tensor, function(x) x[v,v])
  lam.tensor <- lam.tensor[v]
  mu <- mu[v]
  Q <-Q[v,v]

  list(lam.tensor=lam.tensor, mu=mu, Q=Q)
}



#' Akaike Information Criterion (AIC)
#'
#' @param vec a vector of likelihood values
#' @param npar number of parameters
#'
#' @description
#' The Akaike Information Criterion (AIC) is a measure used to compare the relative quality
#' of statistical models for a given dataset. It is calculated as:
#'
#' \deqn{AIC = 2k - 2 \ln(L)}
#'
#' @return AIC scores
#' @export
#'
#' @examples
#'
#' get_aic(2768.948, 4)
get_aic <- function(vec, npar){
  2*npar - 2*vec
}



# lik.const.td <- constrain(lik.td, formulae = formula.td)
# arg.const.td <- argnames(lik.const.td)
# arg.const.td
# init <- starting.point.classe_td(phy, k=3, n.epoch=2)
# starting.point <- init[arg.const.td]
# starting.point




# Creates a starting point for ML inference with episodic ClaSSE
#
# @param phy phylogenetic tree
# @param k character state count
# @param n.epoch number of time epochs
# @param eps Ratio of extinction to speciation rates to be used when choosing a starting set of parameters.See diversitree:::starting.point.classe().
#
# @description
# This function uses diversitree:::starting.point.classe() to create the same
# starting parameter values for each time epoch.
#
# @seealso \code{\link[diversitree]{starting.point.classe}}
#
# @return A named vector containing starting parameter values for each time epoch.
# @export
#
# @examples
#
# file_path <- system.file("extdata", "geosse_tb_tree.rds", package = "divextra")
# phy <- readRDS(file_path)
# starting.point.classe_td(phy, k=3, n.epoch=2)
starting.point.classe_td <- function(phy, k, n.epoch, eps){
  start.classe <- diversitree:::starting.point.classe(phy, k, eps = eps)
  argnames_classe <-  diversitree:::default.argnames.classe(k)
  argnames.td <-  diversitree:::argnames_twopart(argnames_classe, n.epoch)
  out <- rep(start.classe, n.epoch)
  names(out) <- argnames.td
  out
}


#' Creates a starting point for episodic ClaSSE with parameter constrains
#'
#' @param lik.const.td constrained likelihood function
#' @param phy tree
#' @param k number of states
#' @param n.epoch number of epochs
#' @param eps Ratio of extinction to speciation rates to be used when choosing a starting set of parameters.See diversitree:::starting.point.classe().
#'
#' @return A named vector containing starting parameter values for each time epoch.
#' @export
#'
#' @examples
#'
#' file_path <- system.file("extdata", "geosse_tb_tree.rds", package = "divextra")
#' phy <- readRDS(file_path)
#' file_path <- system.file("extdata", "geosse3_td_test.yml", package = "divextra")
#' par.categories.td <- read_yaml_pars_td(file_path)
#' formula.td <- make_constraints_sse_td(par.categories.td)
#' lik.td <-make.classe.td(phy, phy$tip.state, k=3, n.epoch=2, control=list(backend = "gslode"))
#' lik.const.td <- diversitree:::constrain(lik.td, formulae = formula.td)
#' starting.point <- init.pars.classe_td(lik.const.td, phy, k=3, n.epoch=2)
init.pars.classe_td <- function(lik.const.td, phy, k, n.epoch, eps=0.5){
  arg.const.td <- diversitree:::argnames(lik.const.td)
  init <- starting.point.classe_td(phy, k=k, n.epoch=n.epoch, eps=eps)
  starting.point <- init[arg.const.td]
  starting.point
}




#' Group parameters based on the same parameter estimates
#'
#' @param vec A named vector of parameter estimates
#'
#' @return a list
#' @export
#'
pars2groups <- function(vec){
  split(names(vec), vec)
}
