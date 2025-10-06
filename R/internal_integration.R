#' Returns a n*m matrix with the integrals of \eqn{fun_i(s)} along the augmented
#' Turnbull intervals \eqn{I_j}, for \eqn{i=1,\dots,n} and \eqn{j=1,\dots,m}.
#' 
#' Numerical integration is attempted with \code{stats::integrate()} and, on failure,
#' retried with \code{cubature::cubintegrate()}.
#'
#' @param fun Function of the form \code{fun(s, i, j, ...)} returning the
#'   integrand value at scalar \code{s} for subject index \code{i} and interval
#'   index \code{j}. 
#' @param intmap,Kappas computed in gelc.fit 
#' @param relTol Relative tolerance passed to \code{cubature::cubintegrate()} (default \code{1e-5}).
#' @param rel.tol Relative tolerance passed to \code{stats::integrate()} (default \code{relTol}).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{results}}{Numeric \eqn{n \times m} matrix of integrals.}
#'   \item{\code{errmat}}{Numeric \eqn{n \times m} matrix of absolute errors reported by the integrators.}
#'   \item{\code{messmat}}{Character \eqn{n \times m} matrix of status messages (\code{"OK"}, \code{"maxEval reached"}, or \code{"error"}).}
#'   \item{\code{warn}}{Character scalar with aggregated warnings (or \code{NA} if none).}
#' }
#'
#' @seealso \code{\link[stats]{integrate}}, \code{\link[cubature]{cubintegrate}}
#'
#' @importFrom cubature cubintegrate
#' @keywords internal
#' @noRd
integration <- function(fun = function(s, i, ...) dfam(y[i], thetHat[1] + thetHat[2]*s, thetHat[3], famly, size),
                        intmap, Kappas, 
                        relTol = 1e-05, rel.tol = relTol)
{
  
  m <- ncol(Kappas)
  n <- nrow(Kappas)
  
  results <- matrix(0, nrow=n, ncol=m)
  errmat <- matrix(0, nrow=n, ncol=m)
  messmat <- matrix("OK", nrow=n, ncol=m)
  
  # Evaluate the integrand when intervals are exact points
  Ipoint <- abs(intmap[2,] - intmap[1,]) < .Machine$double.eps^0.5
  if(sum(Ipoint)>0) results[,Ipoint] <- sapply(which(Ipoint), function(j) fun(intmap[1,j], 1:n, j))
  
  # Otherwise, perform vectorized numerical approximation of the integral for subjects with kappa=1
  cert <- matrix(TRUE, nrow=n, ncol=m)
  cert[,Ipoint] <- FALSE
  ids <- which(Kappas==1 & cert)
  vect <- cbind((ids - 1) %% n + 1,
                (ids - 1) %/% n + 1)
  if(nrow(vect)>0){
    res <- apply(vect, 1, function(xx) 
      try(stats::integrate(f = function(s) fun(s, xx[1], xx[2]), 
                           lower = intmap[1,xx[2]], 
                           upper = intmap[2,xx[2]],
                           rel.tol = rel.tol),
          silent = TRUE))
    
    # Handle integration failures using cubature
    tryerr <- which(sapply(res, inherits, "try-error"))
    if(length(tryerr)>0){
      for(itry in tryerr){
        xx <- vect[itry,]
        ritry <- cubature::cubintegrate(
          f = function(s) fun(s, xx[1], xx[2]), 
          lower = intmap[1,xx[2]], 
          upper = intmap[2,xx[2]],
          method = "pcubature",
          relTol =  relTol,
          maxEval = 10^12)
        res[[itry]] <- with(ritry, list(value = integral,
                                        abs.error = error,
                                        message = ifelse(returnCode==0, "OK",
                                                  ifelse(returnCode==1, "maxEval reached", "error"))))
      }
    }
    
    # Keep integration results, errors, and messages
    results[ids] <- vapply(res, function(x) x$value, numeric(1))
    errmat[ids] <- vapply(res, function(x) x$abs.error, numeric(1))
    messmat[ids] <- vapply(res, function(x) x$message, character(1))
  }
  
  # Warning messages
  warn = NA
  errs <- ifelse(rowSums(results, na.rm=T)!=0, rowSums(errmat, na.rm=T)/rowSums(results, na.rm=T), 0)
  if(any(errs>0.05 & rowSums(results, na.rm=T)>1e-16)){
    index <- which(errs>1e-3 & rowSums(results, na.rm=T)>1e-16)
    if(is.na(warn)) warn <- deparse(fun)[2]
    for(ii in index)
      warn <- paste0(warn, "\n",
                     "In integrals for subject ", ii, ", for a cum result of ", rowSums(results, na.rm=T)[ii],
                     " the cum error is ", rowSums(errmat, na.rm=T)[ii]," (rel ",errs[ii],")")
  }
  
  if(any(messmat!='OK')){
    index <- which(messmat!='OK')
    row <- (index - 1) %% n + 1
    col <- (index - 1) %/% n + 1
    if(is.na(warn)) warn <- deparse(fun)[2]
    for(ii in 1:length(index))
      warn <- paste0(warn, "\n",
                     "In integral I_", col[ii], " from subj ", row[ii], ": ",
                     messmat[row[ii],col[ii]])
  }
  
  return(mget(c("results","errmat","messmat","warn")))
}