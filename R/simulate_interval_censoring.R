#' Simulate Non-Informative Censoring Intervals Around Latent Variable
#'
#' Generates left (`zl`) and right (`zr`) interval endpoints for each value of a latent variable `z`.
#' The censoring mechanism is non-informative on the true value of the latent variable `z`.
#'
#' @param z Numeric vector. Latent true values.
#' @param zmin Numeric. Lower bound of the support for `z`. Default is 0.
#' @param zmax Optional numeric. Upper bound of the support for `z`. If `NULL`, only `zmin` is enforced.
#' @param rini Optional function. A function that generates the initial gap after `zmin`. If `NULL`, `rfun` is used instead.
#' @param rfun Function. Used to generate random positive gap widths.
#' @param dec Optional integer. Number of decimals to round the generated gaps. If `NULL`, no rounding is applied.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{zl}{Vector of left endpoints of the censoring interval.}
#'   \item{zr}{Vector of right endpoints of the censoring interval.}
#'   \item{gaps}{The full vector of simulated gap widths used to build the intervals.}
#'   \item{warning}{Character string if any simulated interval fails to contain the true `z` value.}
#' }
#'
#' @examples
#' set.seed(1)
#' z <- runif(10, 0, 10)
#' censoring(z, rfun = function(n) rexp(n, 1), dec = 2)
#'
#' @references
#' Gómez, G., Calle, M. L., Oller, R., & Langohr, K. (2009).
#' Tutorial on methods for interval-censored data and their implementation in R.
#' \emph{Statistical Modelling}, 9(4), 259-297. \doi{10.1177/1471082X09009004}
#'
#' @export
censoring <- function(z, zmin = 0, zmax = NULL, rini = NULL, rfun, dec = NULL){
  n <- length(z)
  
  # Check that all z values are within bounds
  if (!is.null(zmax)) {
    if (min(z) < zmin | max(z) > zmax) {
      stop("There are z's outside the support delimited by zmin and zmax")
    }
  } else {
    if (min(z) < zmin) {
      stop("There are z's outside the support delimited by zmin")
    }
  }
  
  # Initial gap from zmin 
  gap_ini <- 0
  if (!is.null(rini)) {
    gap_ini <- rini(n)
    if (!is.null(dec)) gap_ini <- round(gap_ini, dec)
  }
  
  # Function to generate gap widths
  gap_generator <- function(nn) {
    gaps <- rfun(nn)
    if (!is.null(dec)) gaps <- round(gaps, dec)
    return(gaps)
  }
  
  # Accumulate gaps until we exceed each z value
  gaps <- gap_generator(n * ceiling(max(z)))
  mgaps <- matrix(abs(gaps), nrow = n)
  while (any(rowSums(cbind(zmin, gap_ini, mgaps)) < z)) {
    gaps <- c(gaps, gap_generator(n * 5))
    mgaps <- matrix(abs(gaps), nrow = n)
  }
  
  # Cumulative sums to get potential censoring bounds
  U <- t(apply(cbind(zmin, gap_ini, mgaps), MARGIN = 1, FUN = cumsum))
  
  # Determine zl (largest value <= z) and zr (smallest value >= z)
  auxleft <- (z - U > -.Machine$double.eps^0.5)
  auxright <- (U - z > -.Machine$double.eps^0.5)
  
  zl <- apply(auxleft * U, 1, max)
  zr <- apply(auxright * U, 1, function(x) min(x[x > 0]))
  
  # Enforce upper bound if zmax is set
  if (!is.null(zmax)) {
    zr <- pmin(zr, zmax)
  }
  
  # Diagnostic warning
  iwarning <- NULL
  if (any(!(zl - z < .Machine$double.eps^0.5)) & any(!(zr - z > -.Machine$double.eps^0.5))) {
    iwarning <- "Simulated intervals don't contain the true value"
  }
  
  return(list(
    zl = zl,
    zr = zr,
    gaps = gaps,
    warning = iwarning
  ))
}