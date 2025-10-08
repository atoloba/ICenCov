#' Construct an interval-censored covariate for formulas
#'
#' Helper used inside model formulas to bundle left/right endpoints of an
#' interval-censored covariate.
#'
#' @param zl Numeric vector of left endpoints.
#' @param zr Numeric vector of right endpoints.
#' @param name (optional) Name of the variable; defaults to \emph{"z"}.
#'
#' @export
#' @noRd
ic <- function(zl, zr, name = "z"){
  out <- cbind(zl,zr)
  colnames(out) <- paste0(name,c("l","r"))
  return(out)
}
