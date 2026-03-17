#' Normality of residuals for a GLM with an interval-censored covariate
#'
#' Plot-based diagnostic tools to assess normality of interval-censored residuals
#' in generalized linear models with an interval-censored covariate fitted with
#' \code{\link{icglm}}.
#'
#' \code{plot_icres_cdf()} overlays Turnbull's NPMLE-based CDF of the interval-censored
#' residuals with a Normal reference CDF.
#'
#' \code{plot_icres_qq()} produces an interval Q--Q plot by comparing Turnbull-based
#' quantiles of the interval-censored residuals with theoretical Normal quantiles.
#'
#' @param LR A two-column numeric matrix giving interval-censored residuals,
#'   as returned by \code{\link{ires}}.
#' @param hatphi Numeric. Dispersion estimate from \code{\link{icglm}} used to scale
#'   the Normal reference for Pearson and deviance residuals (see Details).
#' @param type Character string specifying the residual type. Options are
#'   \code{"quantile"}, \code{"pearson"}, and \code{"deviance"}.
#'
#' @details
#' Both functions estimate the distribution of the interval-censored residuals using
#' Turnbull's NPMLE via \code{\link[interval]{icfit}}.
#'
#' The Normal reference distribution depends on \code{type}. For Pearson and
#' deviance residuals, the reference is \eqn{N(0, \sqrt{\hat\phi})}, where \eqn{\hat\phi}
#' is provided through \code{hatphi}. For quantile residuals, the reference is \eqn{N(0,1)}.
#'
#' @references
#' Gómez Melis, G., Toloba, A., and Langohr, K. (forthcoming).
#' "Residual Diagnostics for Generalized Linear Models with Interval-Censored Covariates."
#' In: Du, M., Chen, D.-G., Jin, Z., and Sun, J. (eds.),
#' \emph{Next-Gen Lifetime Data Analysis: Emerging Innovations and Applications}.
#' Springer.
#'
#' @importFrom ggplot2 ggplot aes geom_segment geom_abline stat_function labs theme_bw theme
#' @importFrom interval icfit
#' @importFrom stats pnorm qnorm
#'
#' @examples
#' \donttest{
#' ## NOTE: `toy_resnorm` is included with the package as a precomputed toy example.
#' ## The code used to generate it is available in the package source (see `data-raw/toy-icglm-model-diagnostics.R`).
#'
#' data(toy_resnorm)
#'
#' ## This toy dataset was simulated from a Gamma GLM with log link.
#' ## Therefore, the Gamma fit is expected to yield residual diagnostics close to Normal,
#' ## while an inverse Gaussian fit (intentional misspecification) should show clear departures.
#'
#' dat <- toy_resnorm$dat
#' fit_gam <- toy_resnorm$fit_gam
#' fit_ig  <- toy_resnorm$fit_ig
#'
#' ## Plot interval-censored observations
#' plotint(dat)
#'
#' ## Interval-censored residuals
#' X <- model.matrix(y ~ 1, data = dat)
#'
#' res_gam <- ires(
#'   fit_gam$thetaHat, dat$y, dat$zl, dat$zr, X,
#'   Gamma("log"),
#'   type = "all"
#' )
#'
#' res_ig <- ires(
#'   fit_ig$thetaHat, dat$y, dat$zl, dat$zr, X,
#'   inverse.gaussian("log"),
#'   type = "all"
#' )
#'
#' ## Residual diagnostics (quantile residuals)
#' g_gam_ecdf <- plot_icres_cdf(
#'   res_gam$icquant,
#'   hatphi = fit_gam$thetaHat[length(fit_gam$thetaHat)],
#'   type = "quantile"
#' ) + ggplot2::labs(title = NULL, y = "cum. distribution function", x = NULL)
#'
#' g_ig_ecdf <- plot_icres_cdf(
#'   res_ig$icquant,
#'   hatphi = fit_ig$thetaHat[length(fit_ig$thetaHat)],
#'   type = "quantile"
#' ) + ggplot2::labs(title = NULL, y = "cum. distribution function", x = NULL)
#'
#' g_gam_qq <- plot_icres_qq(
#'   res_gam$icquant,
#'   hatphi = fit_gam$thetaHat[length(fit_gam$thetaHat)],
#'   type = "quantile"
#' ) + ggplot2::labs(title = NULL)
#'
#' g_ig_qq <- plot_icres_qq(
#'   res_ig$icquant,
#'   hatphi = fit_ig$thetaHat[length(fit_ig$thetaHat)],
#'   type = "quantile"
#' ) + ggplot2::labs(title = NULL)
#'
#'
#' ## In the output, look for: (i) ECDF close to the dashed Normal curve and
#' ## (ii) Q-Q intervals mostly aligned with the diagonal. Strong deviations suggest lack of fit.
#'
#' ggpubr::ggarrange(g_gam_ecdf, g_gam_qq, g_ig_ecdf, g_ig_qq, ncol = 2, nrow = 2)
#' }
#'
#' @name plot_icres
NULL


#' @rdname plot_icres
#' @export
plot_icres_cdf <- function(LR, hatphi, type) {
  L <- LR[, 1]
  R <- LR[, 2]

  fit_mu <- 0
  fit_sd <- if (type == "quantile") 1 else sqrt(hatphi)

  TBfit <- interval::icfit(L, R, Lin = TRUE, Rin = TRUE)
  pf <- TBfit$pf
  intmap <- TBfit$intmap

  what <- cbind(1, t(intmap), pf, cumsum(pf))
  plot_ic_cdf(list(what), style = "block") +
    stat_function(
      fun = function(x) pnorm(x, mean = fit_mu, sd = fit_sd),
      n = 500, alpha = .5, linetype = 2
    ) +
    labs(title = "Turnbull's NPMLE vs. the normal distribution") +
    theme(legend.position = "none")
}

#' @rdname plot_icres
#' @export
plot_icres_qq <- function(LR, hatphi, type) {
  L <- LR[, 1]
  R <- LR[, 2]

  fit_mu <- 0
  fit_sd <- if (type == "quantile") 1 else sqrt(hatphi)

  TBfit <- interval::icfit(L, R, Lin = TRUE, Rin = TRUE)
  pf <- TBfit$pf
  intmap <- TBfit$intmap
  cdf <- cumsum(pf)

  m <- min(nrow(LR), 100)
  pk <- (seq_len(m) - 0.5) / m
  idx <- sapply(pk, function(pp) which(cdf >= pp)[1])

  ql <- intmap[1, idx]
  qr <- intmap[2, idx]
  zk <- qnorm(pk, mean = fit_mu, sd = fit_sd)

  df <- data.frame(zk, ql, qr)
  dfplot <- subset(df, is.finite(ql) & is.finite(qr))
  dfplot$len <- with(dfplot, qr - ql)

  rnge <- range(c(dfplot$ql, dfplot$qr), na.rm = TRUE)
  eps <- 1e-2 * diff(rnge)
  df2 <- subset(dfplot, len <= eps)
  df2$mid <- (df2$ql + df2$qr) / 2

  ggplot(dfplot, aes(x = zk)) +
    geom_segment(aes(y = ql, yend = qr, xend = zk)) +
    geom_segment(data = df2, aes(y = mid, yend = mid + eps, xend = zk)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    labs(
      x = "Theoretical normal quantiles",
      y = "Turnbull interval quantiles",
      title = "Interval Q-Q plot: Turnbull vs Normal"
    ) +
    theme_bw()
}

utils::globalVariables(c("len", "mid", "zk", "ql", "qr"))
