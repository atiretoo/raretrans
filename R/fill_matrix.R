#' Combine prior and data for transition matrix
#'
#' \code{fill_transition} returns the expected value of the transition
#' matrix combining observed transitions for one time step and a prior
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed transitions.
#' @param P A matrix of the priors for each column. Defaults to uniform.
#' @param priorweight total weight for each column of prior as a percentage of sample size or 1 if negative
#' @param returnType A character vector describing the desired return value. Defaults to "T" the transition matrix.
#'
#' @return The return value depends on parameter returnType.
#' \itemize{
#'   \item A - the summed matrix "filled in" using a dirichlet prior
#'   \item T - just the filled in transition matrix
#'   \item TN - the augmented matrix of fates -- use in calculating the CI or for simulation
#' }
#'
#' @export
#'
#' @examples
#' # Build a simple 3-stage TF list (transition + fertility matrices)
#' T_mat <- matrix(c(0.5, 0.3, 0.0,
#'                   0.2, 0.4, 0.1,
#'                   0.0, 0.1, 0.7), nrow = 3, ncol = 3)
#' F_mat <- matrix(c(0.0, 0.0, 1.5,
#'                   0.0, 0.0, 0.0,
#'                   0.0, 0.0, 0.0), nrow = 3, ncol = 3)
#' TF <- list(T = T_mat, F = F_mat)
#' N  <- c(10, 5, 8)
#'
#' # Default: return filled transition matrix T
#' fill_transitions(TF, N)
#'
#' # Return the full population matrix A = T + F
#' fill_transitions(TF, N, returnType = "A")
#'
#' # Use a prior weight equal to the sample size
#' fill_transitions(TF, N, priorweight = 1)
fill_transitions <- function(TF, N, P = NULL, priorweight = -1, returnType = "T") {
  check_TF(TF)
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]
  if (is.null(P)) {
    # fill in with a uniform prior <- <- <-
    P <- matrix(1 / (order + 1), nrow = order + 1, ncol = order)
  } else {
    if (ncol(P) != order || nrow(P) != (order + 1)) {
      stop("Bad dimensions on P")
    }
  }
  Tfilled <- matrix(NA, nrow = order, ncol = order)
  TN <- matrix(NA, nrow = order + 1, ncol = order)
  for (i in 1:order) {
    observed <- Tmat[, i] * N[i]
    if (priorweight > 0 && N[i] > 0) {
      P[, i] <- P[, i] * priorweight * N[i]
    }
    allfates <- c(observed, N[i] - sum(observed)) + P[, i]
    # missing <- allfates == 0
    # allfates[missing] <- 1
    # allfates[!missing] <- allfates[!missing] + sum(missing)
    Tfilled[, i] <- allfates[1:order] / sum(allfates)
    TN[, i] <- allfates
  }
  if (returnType == "A") {
    return(Tfilled + Fmat)
  } else if (returnType == "T") {
    return(Tfilled)
  } else if (returnType == "TN") {
    return(TN)
  } else {
    stop("Bad returntype in fill_transitions()")
  }
}

#' Combine prior and data for fertility matrix
#'
#' \code{fill_fertility} returns the expected value of the fertility
#' matrix combining observed recruits for one time step and a Gamma prior for each column.
#'
#' Assumes that only one stage reproduces ... needs generalizing.
#'
#' @param TF A list of two matrices, T and F, as ouput by \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed stage distribution.
#' @param alpha A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.
#' @param beta A matrix of the prior parameter for each stage. Impossible stage combinations marked with NA_real_.
#' @param priorweight total weight for each column of prior as a percentage of sample size or 1 if negative
#' @param returnType A character vector describing the desired return value. Defaults to "F" the fertility matrix
#'
#' @return The return value depends on parameter returnType.
#' \itemize{
#'   \item A - the summed matrix "filled in" using a Gamma prior
#'   \item F - just the filled in fertility matrix
#'   \item ab - the posterior parameters alpha and beta as a list.
#' }
#'
#' @export
#'
#' @examples
#' # Build a simple 3-stage TF list (transition + fertility matrices)
#' T_mat <- matrix(c(0.5, 0.3, 0.0,
#'                   0.2, 0.4, 0.1,
#'                   0.0, 0.1, 0.7), nrow = 3, ncol = 3)
#' F_mat <- matrix(c(0.0, 0.0, 1.5,
#'                   0.0, 0.0, 0.0,
#'                   0.0, 0.0, 0.0), nrow = 3, ncol = 3)
#' TF <- list(T = T_mat, F = F_mat)
#' N  <- c(10, 5, 8)
#'
#' # Only adults (stage 3) reproduce; mark non-reproducing entries with NA
#' alpha_mat <- matrix(c(NA, NA, 0.5,
#'                       NA, NA, NA,
#'                       NA, NA, NA), nrow = 3, ncol = 3)
#' beta_mat  <- matrix(c(NA, NA, 1.0,
#'                       NA, NA, NA,
#'                       NA, NA, NA), nrow = 3, ncol = 3)
#'
#' # Default: return filled fertility matrix F
#' fill_fertility(TF, N, alpha = alpha_mat, beta = beta_mat)
#'
#' # Return the full population matrix A = T + F
#' fill_fertility(TF, N, alpha = alpha_mat, beta = beta_mat, returnType = "A")
#'
#' # Return the posterior alpha and beta parameters
#' fill_fertility(TF, N, alpha = alpha_mat, beta = beta_mat, returnType = "ab")
fill_fertility <- function(TF, N, alpha = 0.00001, beta = 0.00001, priorweight = -1, returnType = "F") {
  check_TF(TF)
  Tmat <- TF$T
  Fmat <- TF$F
  order <- dim(Tmat)[1]

  if (length(N) != order || sum(is.na(N)) > 0) {
    stop("N isn't the correct length or has missing values.")
  }

  if ((is.null(dim(alpha)) && length(alpha) != 1) || (is.null(dim(beta)) && length(beta) != 1)) {
    stop("alpha or beta is not a matrix or a single value.")
  }

  if (!(is.numeric(alpha) && is.numeric(beta))) {
    stop("alpha or beta must be numeric matrices or single values.")
  }

  if ((length(alpha) != order^2 || length(beta) != order^2)) {
    warning("length(alpha | beta) != order^2: only using first value of alpha and beta")
    alpha <- matrix(rep(alpha[1], order^2), nrow = order, ncol = order)
    beta <- matrix(rep(beta[1], order^2), nrow = order, ncol = order)
  }

  Ffilled <- matrix(NA, nrow = order, ncol = order)
  babies_next_year <- sweep(Fmat, 2, N, FUN = "*")
  # test reproducing stages with beta, because of divide by zero issues
  reproducing_stages <- apply(beta, 2, function(x) sum(!is.na(x))) > 0
  # matrix multiplication doesn't preserve the column structure

  if ((all(N[reproducing_stages] > 0) | sum(beta, na.rm = TRUE) > 0)) {
    if (priorweight > 0) {
      alpha_post <- sweep(alpha, 2, N, FUN = "*") * priorweight + babies_next_year
      beta_post <- sweep(beta, 2, N, FUN = "*") * priorweight + N
    } else {
      alpha_post <- alpha + babies_next_year
      beta_post <- sweep(beta, 2, N, FUN = "+")
    }
    Ffilled <- alpha_post / beta_post
    Ffilled[is.na(Ffilled)] <- 0
  } else {
    stop("No reproducing stages in N, and no positive values in beta.")
  }

  if (returnType == "A") {
    return(Tmat + Ffilled)
  } else if (returnType == "F") {
    return(Ffilled)
  } else if (returnType == "ab") {
    return(list(alpha = alpha_post, beta = beta_post))
  } else {
    stop("Bad returntype in fill_fertility()")
  }
}

#' Helper functions for generating different priors
#'
#' Extract the number of individuals in each stage from a dataframe
#' of transitions.
#'
#' @param transitions a dataframe of observations of individuals in different stages
#' @param stage the name of the variable with the stage information
#' @param sort a vector of stage names in the desired order. Default is the order of levels in stage.
#'
#' @return a vector of the counts of observations in each level of stage.
#' @export
#'
#' @examples
#' data("L_elto")
#'
#' # Extract one population at one census period
#' onepop <- L_elto[L_elto$POPNUM == 250 & L_elto$year == 5, ]
#' onepop$stage <- factor(onepop$stage, levels = c("p", "j", "a"))
#'
#' # Count individuals per stage in the order p, j, a
#' get_state_vector(onepop, stage = "stage", sort = c("p", "j", "a"))
get_state_vector <- function(transitions, stage = NULL,
                             sort = NULL) {
  if (missing(stage)) {
    stage <- "stage"
  }
  nl <- as.list(1:ncol(transitions))
  names(nl) <- names(transitions)
  stage <- eval(substitute(stage), nl, parent.frame())
  if (is.null(transitions[, stage])) {
    stop("No stage column matching ", stage)
  }
  if (missing(sort)) {
    sort <- levels(transitions[[stage]])
  }
  tf <- table(transitions[[stage]])[sort]
  return(as.vector(tf))
}

check_TF <- function(TF) {

  if(typeof(TF) != "list") {
    stop("TF must be a list of 2 matrices with equal dimensions")
  }
  if(length(TF) != 2) {
    stop("TF must be a list of 2 matrices with equal dimensions")
  }
  if(is.null(TF$T) || is.null(TF$F)) {
    stop("the matrices in TF must be named 'T' and 'F'")
  }
  Tmat <- TF$T
  Fmat <- TF$F

  if(!identical(dim(Tmat), dim(Fmat))) {
    stop("The transition matrix's dimensions don't match the fertility matrix's dimensions")
  }
  if(dim(Tmat)[1] != dim(Tmat)[2]) {
    stop("T and F must be square matrices")
  }

}

