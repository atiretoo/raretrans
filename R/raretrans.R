#' raretrans: Bayesian Priors for Matrix Population Models
#'
#' Provides functions to correct biased transition and fertility estimates in
#' population projection matrices caused by small sample sizes. Small or
#' short-term studies frequently produce structural zeros (biologically possible
#' transitions never observed) and structural ones (transitions estimated at
#' 100% survival, stasis, or mortality that are biologically implausible).
#' Both distort matrix structure and bias estimates of population growth.
#' Functions combine observed transition data with Bayesian prior beliefs to
#' regularise estimates from rare or unobserved events.
#'
#' @importFrom rlang .data
"_PACKAGE"
