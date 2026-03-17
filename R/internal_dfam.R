#' Evaluate density and distribution functions for common GLM families
#'
#' Internal utilities that return pointwise densities/mass functions (\code{dfam})
#' and distribution functions (\code{pfam}) given \code{y}, a GLM family object, a
#' dispersion \code{phi}, and either the mean \code{mu} or linear predictor \code{eta}.
#'
#' @keywords internal
#' @param y Numeric vector of responses.
#' @param eta (optional) Numeric vector: linear predictor on the link scale.
#' @param phi Positive numeric scalar: dispersion parameter. For \code{binomial} and
#'   \code{poisson}, this must be 1.
#' @param famly A \code{stats::family} object.
#' @param size (optional) Parameter for the binomial family; if \code{NULL}, assumes
#'   Bernoulli (\code{size = 1}).
#' @param mu (optional) Numeric vector: mean on the response scale.
#'
#'
#' @details
#' Families supported: \code{gaussian}, \code{Gamma}, \code{inverse.gaussian},
#' \code{poisson}, \code{binomial}. For inverse Gaussian, \pkg{statmod} is used.
#'
#' @name internal_dfam
#' @keywords internal
NULL


#' @rdname internal_dfam
#' @importFrom stats dgamma dnorm dpois dbinom
#' @importFrom statmod dinvgauss
dfam <- function(y, eta = NULL, phi, famly, size = NULL, mu = NULL) {
  if (is.null(mu)) mu <- famly$linkinv(eta)

  if (famly$family == "Gamma") {
    dgamma(y, shape = 1/phi, scale = mu * phi)

  } else if (famly$family == "gaussian") {
    dnorm(y, mean = mu, sd = sqrt(phi))

  } else if (famly$family == "inverse.gaussian") {
    statmod::dinvgauss(y, mean = mu, dispersion = phi)

  } else if (famly$family == "poisson") {
    if (phi != 1) warning("Dispersion must be 1")
    dpois(y, lambda = mu)

  } else if (famly$family == "binomial") {
    if (phi != 1) warning("Dispersion must be 1")
    sze <- if (is.null(size)) 1 else size
    dbinom(y, size = sze, prob = mu / sze)

  } else stop("Unsupported family in dfam()")
}

#' @rdname internal_dfam
#' @importFrom stats pgamma pnorm ppois pbinom
#' @importFrom statmod pinvgauss
pfam <- function(y, eta = NULL, phi, famly, size = NULL, mu = NULL) {
  if (is.null(mu)) mu <- famly$linkinv(eta)

  if (famly$family == "Gamma") {
    pgamma(y, shape = 1/phi, scale = mu * phi)

  } else if (famly$family == "gaussian") {
    pnorm(y, mean = mu, sd = sqrt(phi))

  } else if (famly$family == "inverse.gaussian") {
    statmod::pinvgauss(y, mean = mu, dispersion = phi)

  } else if (famly$family == "poisson") {
    if (phi != 1) warning("Dispersion must be 1")
    ppois(y, lambda = mu)

  } else if (famly$family == "binomial") {
    if (phi != 1) warning("Dispersion must be 1")
    sze <- if (is.null(size)) 1 else size
    pbinom(y, size = sze, prob = mu / sze)

  } else stop("Unsupported family in pfam()")
}
