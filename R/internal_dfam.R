#' Evaluate densities for common GLM families
#'
#' Internal utility that returns pointwise densities (or mass functions for
#' discrete families) given `y`, a GLM family object, a dispersion `phi`, and
#' either the mean `mu` or the linear predictor `eta`.
#'
#' @keywords internal
#' @param y Numeric vector of responses.
#' @param eta (optional) Numeric vector: linear predictor on the link scale.
#' @param phi Positive numeric scalar: dispersion parameter. For
#'   \code{binomial} and \code{poisson}, this must be 1.
#' @param famly A \code{stats::family} object.
#' @param size (optional) Parameter for the binomial family; if `NULL`, assumes Bernoulli (logistic, `size = 1`).
#' @param mu (optional) Numeric vector: mean on the response scale.
#'
#' @return Numeric vector of densities (or probabilities for discrete
#'   families), same length as `y`.
#'
#' @details
#' \itemize{
#'   \item \strong{Gaussian}: \code{dnorm(y, mean = mu, sd = sqrt(phi))}.
#'   \item \strong{Gamma}: \code{dgamma(y, shape = 1/phi, scale = mu * phi)}.
#'   \item \strong{Inverse Gaussian}: \code{statmod::dinvgauss(y, mean = mu, dispersion = phi)}.
#'   \item \strong{Binomial}: \code{dbinom(y, size = size, prob = mu / size)} (requires \code{phi = 1}).
#'   \item \strong{Poisson}: \code{dpois(y, lambda = mu)} (requires \code{phi = 1}).
#' }
#'
#' @importFrom statmod dinvgauss
#' @keywords internal
#' @noRd
dfam <- function(y, eta=NULL, phi, famly, size = NULL, mu=NULL) {
  if(is.null(mu)) mu <- famly$linkinv(eta)
  if(famly$family=="Gamma") {
    dgamma(y, shape = 1/phi, scale = mu*phi)
    
  } else if(famly$family == "binomial"){
    if(phi != 1) warning("Dispersion must be 1")
    if(is.null(size)) sze = 1 else sze = size
    dbinom(y, size = sze, prob = mu/sze)
    
  } else if(famly$family == "gaussian"){
    dnorm(y, mean = mu, sd = sqrt(phi))
    
  } else if(famly$family == "inverse.gaussian"){
    statmod::dinvgauss(y, mean = mu, dispersion = phi)
    
  } else if(famly$family == "poisson"){
    if(phi != 1) warning("Dispersion must be 1")
    dpois(y, lambda = mu)
  }
}