#' Calculate beta credible intervals for all transition matrix entries
#'
#' Computes the marginal posterior beta credible intervals for every entry
#' in the transition matrix, including the probability of dying. The marginal
#' posterior distribution of each transition probability follows a beta
#' distribution derived from the Dirichlet-multinomial model, using the
#' augmented fate matrix (TN) returned by \code{\link{fill_transitions}}.
#'
#' @param TF A list of two matrices, T and F, as output by
#'   \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed stage distribution at the start of the
#'   transition period.
#' @param P A matrix of priors for each column. Defaults to uniform.
#' @param priorweight Total weight for each column of prior as a percentage
#'   of sample size, or 1 if negative. Defaults to -1 (uninformative).
#' @param ci Credible interval width as a probability between 0 and 1.
#'   Defaults to 0.95 (95\% credible interval).
#' @param stage_names Optional character vector of stage names in the same
#'   order as the columns of the transition matrix. If \code{NULL}, names
#'   are taken from \code{colnames(TF$T)}, or generic labels
#'   \code{"Stage 1"}, \code{"Stage 2"}, etc. are used.
#'
#' @return A data frame with one row per fate per stage (including the dead
#'   fate) and the following columns:
#' \describe{
#'   \item{from_stage}{Character. The source stage (column of the matrix).}
#'   \item{to_stage}{Character. The destination stage, including \code{"dead"}.}
#'   \item{mean}{Numeric. Posterior mean transition probability.}
#'   \item{lower}{Numeric. Lower bound of the credible interval.}
#'   \item{upper}{Numeric. Upper bound of the credible interval.}
#' }
#' @export
#'
#' @seealso \code{\link{plot_transition_CrI}} for visualising the output.
#'
#' @examples
#' T_mat <- matrix(c(0.5, 0.3, 0.0,
#'                   0.2, 0.4, 0.1,
#'                   0.0, 0.1, 0.7), nrow = 3, ncol = 3)
#' F_mat <- matrix(c(0.0, 0.0, 1.5,
#'                   0.0, 0.0, 0.0,
#'                   0.0, 0.0, 0.0), nrow = 3, ncol = 3)
#' TF <- list(T = T_mat, F = F_mat)
#' N  <- c(10, 5, 8)
#'
#' # Default 95% credible intervals
#' cri <- transition_CrI(TF, N, stage_names = c("plantula", "juvenile", "adult"))
#' cri
#'
#' # 90% credible intervals
#' transition_CrI(TF, N, stage_names = c("plantula", "juvenile", "adult"), ci = 0.90)
transition_CrI <- function(TF, N, P = NULL, priorweight = -1, ci = 0.95,
                           stage_names = NULL) {
  check_TF(TF)

  if (!is.numeric(ci) || length(ci) != 1 || ci <= 0 || ci >= 1) {
    stop("ci must be a single numeric value strictly between 0 and 1.")
  }

  order <- dim(TF$T)[1]

  if (is.null(stage_names)) {
    if (!is.null(colnames(TF$T))) {
      stage_names <- colnames(TF$T)
    } else {
      stage_names <- paste0("Stage ", seq_len(order))
    }
  } else {
    if (length(stage_names) != order) {
      stop("stage_names must have the same length as the number of stages.")
    }
  }

  to_names    <- c(stage_names, "dead")
  alpha_lower <- (1 - ci) / 2
  alpha_upper <- 1 - alpha_lower

  TN <- fill_transitions(TF, N, P = P, priorweight = priorweight,
                         returnType = "TN")

  results <- vector("list", order)
  for (j in seq_len(order)) {
    a   <- TN[, j]
    b   <- sum(TN[, j]) - a
    p   <- a / (a + b)
    lcl <- stats::qbeta(alpha_lower, a, b)
    ucl <- stats::qbeta(alpha_upper, a, b)

    results[[j]] <- data.frame(
      from_stage = stage_names[j],
      to_stage   = to_names,
      mean       = p,
      lower      = lcl,
      upper      = ucl,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}


#' Plot beta credible intervals for transition matrix entries
#'
#' Creates a \code{ggplot2} visualisation of the posterior mean transition
#' probabilities and their credible intervals for each stage, including the
#' probability of dying. Each panel shows the fate distribution from one
#' source stage.
#'
#' @param cri A data frame as returned by \code{\link{transition_CrI}}.
#' @param include_dead Logical. Whether to include the dead fate in the plot.
#'   Defaults to \code{TRUE}.
#' @param title Character. Plot title. Defaults to
#'   \code{"Posterior transition probabilities with credible intervals"}.
#'
#' @return A \code{ggplot} object.
#' @export
#'
#' @seealso \code{\link{transition_CrI}} for computing the credible intervals.
#'
#' @examples
#' T_mat <- matrix(c(0.5, 0.3, 0.0,
#'                   0.2, 0.4, 0.1,
#'                   0.0, 0.1, 0.7), nrow = 3, ncol = 3)
#' F_mat <- matrix(c(0.0, 0.0, 1.5,
#'                   0.0, 0.0, 0.0,
#'                   0.0, 0.0, 0.0), nrow = 3, ncol = 3)
#' TF <- list(T = T_mat, F = F_mat)
#' N  <- c(10, 5, 8)
#'
#' cri <- transition_CrI(TF, N, stage_names = c("plantula", "juvenile", "adult"))
#'
#' # Include dead fate (default)
#' plot_transition_CrI(cri)
#'
#' # Exclude dead fate
#' plot_transition_CrI(cri, include_dead = FALSE)
plot_transition_CrI <- function(cri, include_dead = TRUE,
                                title = "Posterior transition probabilities with credible intervals") {

  if (!is.data.frame(cri) ||
      !all(c("from_stage", "to_stage", "mean", "lower", "upper") %in% names(cri))) {
    stop("cri must be a data frame as returned by transition_CrI().")
  }

  if (!include_dead) {
    cri <- cri[cri$to_stage != "dead", ]
  }

  cri$from_stage <- factor(cri$from_stage, levels = unique(cri$from_stage))
  cri$to_stage   <- factor(cri$to_stage,   levels = unique(cri$to_stage))

  ggplot2::ggplot(cri,
                  ggplot2::aes(x    = .data$to_stage,
                               y    = .data$mean,
                               ymin = .data$lower,
                               ymax = .data$upper)) +
    ggplot2::geom_pointrange() +
    ggplot2::facet_wrap(~from_stage, labeller = ggplot2::label_both) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x     = "Destination stage",
      y     = "Transition probability",
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}


#' Plot posterior beta density curves for all transition matrix entries
#'
#' Creates a \code{ggplot2} visualisation arranged as an \eqn{n \times n}
#' matrix of density plots, mirroring the structure of the population
#' projection matrix. Each panel shows the full marginal posterior beta
#' distribution for one transition probability, with the credible interval
#' region shaded. Columns correspond to source stages (from) and rows to
#' destination stages (to), including the dead fate as the bottom row.
#'
#' @param TF A list of two matrices, T and F, as output by
#'   \code{\link[popbio]{projection.matrix}}.
#' @param N A vector of observed stage distribution at the start of the
#'   transition period.
#' @param P A matrix of priors for each column. Defaults to uniform.
#' @param priorweight Total weight for each column of prior as a percentage
#'   of sample size, or 1 if negative. Defaults to -1 (uninformative).
#' @param ci Credible interval width as a probability between 0 and 1.
#'   The shaded region in each panel covers this interval. Defaults to 0.95.
#' @param stage_names Optional character vector of stage names in the same
#'   order as the columns of the transition matrix. If \code{NULL}, names
#'   are taken from \code{colnames(TF$T)}, or generic labels are used.
#' @param include_dead Logical. Whether to include the dead fate as the bottom
#'   row of the matrix plot. Defaults to \code{TRUE}.
#' @param title Character. Plot title.
#'
#' @return A \code{ggplot} object arranged as an \eqn{n \times n} grid (or
#'   \eqn{(n+1) \times n} when \code{include_dead = TRUE}), with source stages
#'   as columns and destination stages as rows.
#' @export
#'
#' @seealso \code{\link{transition_CrI}}, \code{\link{plot_transition_CrI}}.
#'
#' @examples
#' T_mat <- matrix(c(0.5, 0.3, 0.0,
#'                   0.2, 0.4, 0.1,
#'                   0.0, 0.1, 0.7), nrow = 3, ncol = 3)
#' F_mat <- matrix(c(0.0, 0.0, 1.5,
#'                   0.0, 0.0, 0.0,
#'                   0.0, 0.0, 0.0), nrow = 3, ncol = 3)
#' TF <- list(T = T_mat, F = F_mat)
#' N  <- c(10, 5, 8)
#'
#' # Include dead fate as bottom row (default)
#' plot_transition_density(TF, N,
#'                         stage_names = c("plantula", "juvenile", "adult"))
#'
#' # Transitions only, no dead row
#' plot_transition_density(TF, N,
#'                         stage_names = c("plantula", "juvenile", "adult"),
#'                         include_dead = FALSE)
plot_transition_density <- function(TF, N, P = NULL, priorweight = -1,
                                    ci = 0.95, stage_names = NULL,
                                    include_dead = TRUE,
                                    title = "Posterior transition probability densities") {
  check_TF(TF)

  if (!is.numeric(ci) || length(ci) != 1 || ci <= 0 || ci >= 1) {
    stop("ci must be a single numeric value strictly between 0 and 1.")
  }

  order <- dim(TF$T)[1]

  if (is.null(stage_names)) {
    if (!is.null(colnames(TF$T))) {
      stage_names <- colnames(TF$T)
    } else {
      stage_names <- paste0("Stage ", seq_len(order))
    }
  } else {
    if (length(stage_names) != order) {
      stop("stage_names must have the same length as the number of stages.")
    }
  }

  to_names    <- c(stage_names, "dead")
  alpha_lower <- (1 - ci) / 2
  alpha_upper <- 1 - alpha_lower

  TN <- fill_transitions(TF, N, P = P, priorweight = priorweight,
                         returnType = "TN")

  x_seq   <- seq(0, 1, length.out = 300)
  results <- vector("list", order * (order + 1))
  idx     <- 1L

  for (j in seq_len(order)) {
    col_sum <- sum(TN[, j])
    for (i in seq_len(order + 1)) {
      a   <- TN[i, j]
      b   <- col_sum - a
      lcl <- stats::qbeta(alpha_lower, a, b)
      ucl <- stats::qbeta(alpha_upper, a, b)

      dens              <- stats::dbeta(x_seq, a, b)
      dens[!is.finite(dens)] <- NA

      ci_dens           <- dens
      ci_dens[x_seq < lcl | x_seq > ucl] <- 0

      results[[idx]] <- data.frame(
        from_stage = stage_names[j],
        to_stage   = to_names[i],
        x          = x_seq,
        density    = dens,
        ci_density = ci_dens,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  df            <- do.call(rbind, results)

  if (!include_dead) {
    df       <- df[df$to_stage != "dead", ]
    to_names <- to_names[to_names != "dead"]
  }

  df$from_stage <- factor(df$from_stage, levels = stage_names)
  df$to_stage   <- factor(df$to_stage,   levels = to_names)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = .data$ci_density),
      fill = "steelblue", alpha = 0.4, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$density),
      na.rm = TRUE
    ) +
    ggplot2::facet_grid(to_stage ~ from_stage,
                        labeller = ggplot2::labeller(
                          from_stage = ggplot2::label_value,
                          to_stage   = ggplot2::label_value
                        )) +
    ggplot2::labs(
      x     = "Transition probability",
      y     = "Density",
      title = title
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.text  = ggplot2::element_text(size = 8, face = "bold"),
      axis.text   = ggplot2::element_text(size = 6),
      axis.title  = ggplot2::element_text(size = 9),
      panel.spacing = ggplot2::unit(0.3, "lines")
    )
}
