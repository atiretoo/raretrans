#' Bootstrap observed census transitions
#'
#' Calculate bootstrap distributions of population growth rates (lambda), stage vectors, and projection matrix elements by randomly sampling with replacement from a stage-fate data frame of observed transitions
#'
#' @param transitions a stage-fate data frame with stage or age class in the current census, fate in the subsequent census, and one or more fertility columns
#' @param iterations Number of bootstrap iterations
#' @param by.stage.counts Resample transitions with equal probability (default) or by subsets of initial stage counts
#' @param ... additional options passed to \code{\link[raretrans]{projection_matrix}}
#'
#' @return A list with 3 items
#'
#' lambda	A vector containing bootstrap values for lambda
#'
#' matrix	A matrix containing bootstrap transtion matrices with one projection matrix per row.
#'
#' vector	A matrix containing bootstrap stage vectors with one stage vector per row.
#'
#' @export
#'
#' @source This is a modified version of \code{\link[popbio]{boot.transitions}}.

boot_transitions <- function (transitions, iterations, by.stage.counts = FALSE,
          ...)
{
  t <- iterations
  mat <- vector("list", t)
  vec <- vector("list", t)
  lam <- numeric(t)
  for (i in 1:t) {
    if (by.stage.counts) {
      boot <- do.call(rbind, lapply(split(transitions,
                                          transitions$stage, drop = TRUE), function(x) x[sample(nrow(x),
                                                                                                replace = TRUE), ]))
    }
    else {
      boot <- transitions[sample(nrow(transitions), replace = TRUE),
      ]
    }
    A <- projection_matrix(boot, ...)
    vec[[i]] <- table(boot$stage)
    mat[[i]] <- as.vector(A)
    lam[i] <- lambda(A)
  }
  n <- dim(A)[1]
  boot.stage <- list(lambda = lam, matrix = matrix(unlist(mat),
                                                   byrow = TRUE, nrow = t, dimnames = list(1:t, paste("a",
                                                                                                      1:n, rep(1:n, each = n), sep = ""))), vector = matrix(unlist(vec),                                                                                                                                                            byrow = TRUE, nrow = t, dimnames = list(1:t, names(vec[[1]]))))
  boot.stage
}
