
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



#Nstates = 2
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
# # Test
# Nstates = 3
# argsHiClaSSE2 <- argnames_HiClaSSE(Nstates)
# argsHiClaSSE2$pars
# argsHiClaSSE2$arrays
# pars <- c(1:27)
# names(pars) <- argsHiClaSSE2$pars
# pars
# pars_to_arrays(pars, Nstates)
