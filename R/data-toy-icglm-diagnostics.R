#' Toy examples for model diagnostics in a GLM with an interval-censored covariate
#'
#' These objects are simulated toy examples shipped with the package to
#' illustrate diagnostic checks for \code{icglm} models without requiring long
#' model fitting in \code{@examples}.
#'
#' The toy examples correspond to those in:
#' Gómez Melis, G., Toloba, A., and Langohr, K. (forthcoming).
#' "Residual Diagnostics for Generalized Linear Models with Interval-Censored Covariates."
#' In: Du, M., Chen, D.-G., Jin, Z., and Sun, J. (eds.),
#' \emph{Next-Gen Lifetime Data Analysis: Emerging Innovations and Applications}.
#' Springer.
#'
#' @details
#' Each object is a named list that contains a dataset plus pre-fitted \code{icglm} models.
#'
#' @name toy_icglm_diagnostics
NULL

#' @rdname toy_icglm_diagnostics
"toy_resnorm"

#' @rdname toy_icglm_diagnostics
"toy_linearity"

#' @rdname toy_icglm_diagnostics
"toy_link"
