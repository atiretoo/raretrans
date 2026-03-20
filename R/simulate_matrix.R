#' Dirichlet distributed random numbers
#'
#' Generate n random vectors distributed according to a Dirichlet distribution.
#' Each row of the returned matrix is a random vector that sums to 1.
#'
#' @param n The number of random vectors to generate
#' @param alpha A vector of parameters
#'
#' @return The function returns a matrix with n rows and \code{length(alpha)} columns
#' @export
#'
#' @examples
#' # Generate 5 random probability vectors from a Dirichlet(1, 2, 3) distribution
#' set.seed(42)
#' rdirichlet(5, c(1, 2, 3))
#'
#' @source copied from package \code{MCMCpack} to avoid a dependency. That code
#' was taken from Greg's Miscellaneous Functions (gregmisc). His code was based
#' on code posted by Ben Bolker to R-News on 15 Dec 2000.
#'
rdirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(stats::rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  return(x / as.vector(sm))
}

#' Simulate population projection matrices
#'
#' \code{sim_transition} generates a list of simulated population projection matrices from the provided parameters
#' and prior distributions.
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed stages at start of transition.
#' @param P A matrix of the priors for each column. Defaults to uniform.
#' @param alpha A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.
#' @param beta A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.
#' @param priorweight total weight for each column of prior as a percentage of sample size or 1 if negative
#' @param samples The number of matrices to return.
#'
#' @return Always returns a list.
#' @export
#'
#' @examples
#' # Build a simple 3-stage TF list
#' T_mat <- matrix(c(0.5, 0.3, 0.0,
#'                   0.2, 0.4, 0.1,
#'                   0.0, 0.1, 0.7), nrow = 3, ncol = 3)
#' F_mat <- matrix(c(0.0, 0.0, 1.5,
#'                   0.0, 0.0, 0.0,
#'                   0.0, 0.0, 0.0), nrow = 3, ncol = 3)
#' TF <- list(T = T_mat, F = F_mat)
#' N  <- c(10, 5, 8)
#'
#' # Simulate 10 population projection matrices using default uninformative priors
#' set.seed(42)
#' mats <- sim_transitions(TF, N, samples = 10)
#' length(mats)   # 10 matrices
#' mats[[1]]      # first simulated matrix
sim_transitions <- function(TF, N, P = NULL, alpha = 0.00001, beta = 0.00001, priorweight = -1, samples = 1) {
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]
  if (is.null(P)) {
    # fill in with a uniform prior <- <- <-
    P <- matrix(1 / (order + 1), nrow = order + 1, ncol = order)
  } else {
    if (ncol(P) != order | nrow(P) != (order + 1)) {
      stop("Bad dimensions on P")
    }
  }
  TN <- fill_transitions(TF, N, P, priorweight, returnType = "TN")
  # alpha and beta are checked in fill_fertility
  ab_post <- fill_fertility(TF, N, alpha = alpha, beta = beta, priorweight = priorweight, returnType = "ab")
  alpha <- ab_post$alpha # these will now be square matrices
  beta <- ab_post$beta # these will now be square matrices

  Amats <- list()
  for(s in 1:samples){
    T_ <- matrix(0, nrow=order, ncol=order)
    F_ <- matrix(0, nrow=order, ncol=order)

    for(j in 1:order){ # looping over the column of the matrix
      T_[,j] <- rdirichlet(1, TN[,j])[1:order]
      for (i in 1:order){ # loop over rows as well to generate fertilities
        if (!is.na(alpha[i,j])){
          F_[i,j] <- stats::rgamma(1, shape = alpha[i,j], rate = beta[i,j])
        }
      }

    }
    Amats[[s]] <- T_ + F_
  }


  return(Amats)
}
