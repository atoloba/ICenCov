#' Helpers for \code{icglm} objects
#'
#' @description
#' Minimal S3 methods to work with model fits of class \code{icglm},
#' analogous to \code{glm} methods

#'
#' @seealso
#' \code{\link[stats]{coef.glm}}, \code{\link[stats]{vcov.glm}},
#' \code{\link[stats]{summary.glm}}
#'
#' @examples
#' fit <- icglm(RNA ~ age + ic(zl, zr, "time"),
#'              family = Gamma("log"), data = actg359,
#'              Lin = FALSE, Rin = TRUE)
#' summary(fit)
#' coef(fit)
#' vcov(fit)
#'
#' @name icglm-methods
#' @aliases
#' coef.icglm
#' dispersion.icglm
#' vcov.icglm
#' summary.icglm
#' print.summary.icglm
#' @family icglm
NULL

#' @rdname icglm-methods
#' @export
#' @method coef icglm
coef.icglm <- function(object, ...) {
  object$thetaHat
}

#' @rdname icglm-methods
#' @export
dispersion.icglm <- function(object, ...) {
  that <- object$thetaHat
  if (!is.null(names(that)) && "disp" %in% names(that)) unname(that[["disp"]])
  else NA
}

#' @rdname icglm-methods
#' @export
#' @method vcov icglm
vcov.icglm <- function(object, ...) {
  object$vcov
}

#' @rdname icglm-methods
#' @export
#' @method summary icglm
summary.icglm <- function(object, ...) {
  
  co <- coef(object)
  if (!is.null(names(co)) && "disp" %in% names(co)) {
    co <- co[setdiff(names(co), "disp")]
  }
  
  V  <- vcov(object)
  rn <- rownames(V)
  if (!is.null(rn) && "disp" %in% rn) {
    keep <- setdiff(rn, "disp")
    V <- V[keep, keep, drop = FALSE]
  }
  
  se <- sqrt(diag(V))
  z  <- co / se
  p  <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  coef.tab <- cbind(
    Estimate     = co,
    `Std. Error` = se,
    `z value`    = z,
    `Pr(>|z|)`   = p
  )
  
  out <- list(
    call         = object$call,
    family       = object$family,
    coefficients = coef.tab,
    cov          = V,
    dispersion   = dispersion.icglm(object)
  )
  class(out) <- "summary.icglm"
  out
}

#' @rdname icglm-methods
#' @export
#' @method print summary.icglm
print.summary.icglm <- function(x,
                                digits = max(3, getOption("digits") - 3),
                                signif.stars = getOption("show.signif.stars"),
                                ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nFamily: ", x$family$family, "(", x$family$link, ")\n", sep = "")
  cat("\nDispersion:", format(x$dispersion, digits = digits), "\n")
  
  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars)
  invisible(x)
}

