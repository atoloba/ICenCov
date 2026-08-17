#' Asset for model diagnostics of GLMs with an interval-censored covariate
#'
#' Supports model diagnostics when interval uncertainty in an interval-censored covariate
#' propagates into diagnostic plot coordinates for \code{\link{icglm}} models, by mapping each
#' interval through user-supplied functions \code{funx} and \code{funy}.
#'
#' @param LR A two-column numeric matrix giving interval-censored values.
#' @param funx Function applied to grid points within each interval to produce the x-axis
#'   values. Defaults to the identity function.
#' @param funy Function applied to grid points within each interval to produce the y-axis
#'   values.
#' @param lseq Integer. Number of grid points used within each interval. Larger values
#'   produce smoother curves but increase computation time.
#' @param eps Optional non-negative numeric value. If both the horizontal and vertical
#'   ranges of a mapped interval are smaller than or equal to \code{eps}, the interval
#'   is represented by a point instead of a path. Defaults to \code{NULL}, in which case
#'   all mapped intervals are represented by paths.
#' @param funx_args Named list of additional arguments passed to \code{funx}.
#' @param funy_args Named list of additional arguments passed to \code{funy}.
#' @param path_args Named list of additional arguments passed to
#'   \code{\link[ggplot2]{geom_path}}, such as \code{color},
#'   \code{linewidth}, \code{linetype}, and \code{alpha}. Defaults to
#'   \code{list(alpha = 0.4)}.
#' @param point_args Named list of additional arguments passed to
#'   \code{\link[ggplot2]{geom_point}} when \code{eps} is not \code{NULL},
#'   such as \code{color}, \code{fill}, \code{size}, \code{shape},
#'   \code{alpha}, and \code{stroke}. Defaults to
#'   \code{list(alpha = 0.7)}.
#'
#' @details
#' For convenience, the returned \code{ggplot2} object has an attribute \code{"plotLR"},
#' a data frame summarizing the mapped ranges per observation:
#' \code{xl, xr} for the x-range and \code{yl, yr} for the y-range (when applicable).
#'
#'
#' @references
#' Gómez Melis, G., Toloba, A., and Langohr, K. (2026).
#' "Residual Diagnostics for Generalized Linear Models with Interval-Censored Covariates."
#' In: Du, M., Chen, D.-G., Jin, Z., and Sun, J. (eds.),
#' \emph{Next-Gen Lifetime Data Analysis: Emerging Innovations and Applications}.
#' Springer.
#'
#'
#' @importFrom ggplot2 ggplot aes geom_path geom_point theme_bw
#'
#' @examples
#' \donttest{
#' ## NOTE: `toy_linearity` and `toy_link` are included with the package as precomputed toy examples.
#' ## The code used to generate them is available in the package source
#' ## (see `data-raw/toy-icglm-model-diagnostics.R` in the development repository).
#'
#' ## ------------------------------------------------------------
#' ## Linearity example: partial residuals for covariates
#' ## ------------------------------------------------------------
#' data(toy_linearity)
#'
#' ## This toy dataset was simulated so that the true mean structure involves x^2 and log(z).
#' ## The transformed fit is therefore expected to behave better than the naive linear fit.
#'
#' dat <- toy_linearity$dat
#' fit_linear <- toy_linearity$fit_linear
#' fit_trans  <- toy_linearity$fit_trans
#'
#' ## Plot interval-censored observations
#' plotint(dat)
#'
#' ## --- Transformed model: partial residuals vs x^2 and log(z) ---
#' mod <- fit_trans
#' y <- dat$y
#' x2 <- dat$x2
#' zl2 <- dat$zl2
#' zr2 <- dat$zr2
#' theta <- mod$thetaHat
#' famly <- mod$family
#'
#' f.eta  <- function(s) theta[1] + theta[2] * x2 + theta[3] * s
#' f.mu   <- function(s) famly$linkinv(f.eta(s))
#' f.zeta <- function(s) famly$linkfun(f.mu(s)) + (y - f.mu(s)) / famly$mu.eta(f.eta(s))
#'
#' ## Partial residual component for x^2
#' f.res_x2 <- function(s) (f.zeta(s) - f.eta(s)) + theta[2] * x2
#' g_x2 <- plot_ic_mapping(cbind(zl2, zr2),
#'   funx = function(s) x2, funy = f.res_x2
#' ) + ggplot2::labs(x = expression(x^2), y = "partial residuals", title = NULL)
#'
#' plr <- attr(g_x2, "plotLR")
#' yl <- stats::quantile(c(plr$yl, plr$yr), 0, na.rm = TRUE)
#' yu <- stats::quantile(c(plr$yl, plr$yr), 0.97, na.rm = TRUE)
#' g_x2 <- g_x2 + ggplot2::coord_cartesian(ylim = c(yl, yu))
#'
#' ## Partial residual component for log(z)
#' f.res_logz <- function(s) (f.zeta(s) - f.eta(s)) + theta[3] * s
#' g_logz <- plot_ic_mapping(cbind(zl2, zr2),
#'   funx = function(s) s, funy = f.res_logz, lseq = 50
#' ) + ggplot2::labs(x = expression(log(z)), y = "partial residuals", title = NULL)
#'
#' plr <- attr(g_logz, "plotLR")
#' yl <- stats::quantile(c(plr$yl, plr$yr), 0, na.rm = TRUE)
#' yu <- stats::quantile(c(plr$yl, plr$yr), 0.97, na.rm = TRUE)
#'
#' d <- ggplot2::ggplot_build(g_logz)$data[[1]]
#' d_in <- d[is.finite(d$x) & is.finite(d$y) & d$y >= yl & d$y <= yu, ]
#' xl <- min(d_in$x, na.rm = TRUE)
#' xu <- max(d_in$x, na.rm = TRUE)
#' g_logz <- g_logz + ggplot2::coord_cartesian(ylim = c(yl, yu), xlim = c(xl, xu))
#'
#' ## --- Naive model: partial residuals vs x and z (raw scale) ---
#' mod <- fit_linear
#' y <- dat$y
#' x <- dat$x
#' zl <- dat$zl
#' zr <- dat$zr
#' theta <- mod$thetaHat
#' famly <- mod$family
#'
#' f.eta  <- function(s) theta[1] + theta[2] * x + theta[3] * s
#' f.mu   <- function(s) famly$linkinv(f.eta(s))
#' f.zeta <- function(s) famly$linkfun(f.mu(s)) + (y - f.mu(s)) / famly$mu.eta(f.eta(s))
#'
#' f.res_x <- function(s) (f.zeta(s) - f.eta(s)) + theta[2] * x
#' g_x <- plot_ic_mapping(cbind(zl, zr),
#'   funx = function(s) x, funy = f.res_x
#' ) + ggplot2::labs(x = "x", y = "partial residuals", title = NULL)
#'
#' f.res_z <- function(s) (f.zeta(s) - f.eta(s)) + theta[3] * s
#' g_z <- plot_ic_mapping(cbind(zl, zr),
#'   funx = function(s) s, funy = f.res_z
#' ) + ggplot2::labs(x = "z", y = "partial residuals", title = NULL)
#'
#' ggpubr::ggarrange(plotlist = list(g_x2, g_logz, g_x, g_z), nrow = 2, ncol = 2)
#'
#' ## ------------------------------------------------------------
#' ## Link-function example: working response vs linear predictor
#' ## ------------------------------------------------------------
#' data(toy_link)
#'
#' ## This toy dataset was simulated from a Gamma GLM with log link.
#' ## Comparing the identity-link fit vs the log-link fit illustrates how the
#' ## working-response mapping changes under link misspecification.
#'
#' dat <- toy_link$dat
#' fit_id  <- toy_link$fit_id
#' fit_log <- toy_link$fit_log
#'
#' y <- dat$y
#' zl <- dat$zl
#' zr <- dat$zr
#'
#' ## Identity-link fit
#' mod <- fit_id
#' theta <- mod$thetaHat
#' famly <- mod$family
#' f.eta  <- function(s) theta[1] + theta[2] * s
#' f.mu   <- function(s) famly$linkinv(f.eta(s))
#' f.zeta <- function(s) famly$linkfun(f.mu(s)) + (y - f.mu(s)) / famly$mu.eta(f.eta(s))
#'
#' g0 <- plot_ic_mapping(cbind(zl, zr), funx = f.eta, funy = f.zeta) +
#'   ggplot2::geom_abline(linetype = 2, alpha = .6) +
#'   ggplot2::labs(x = "linear predictor", y = "working response", title = NULL)
#' plr <- attr(g0, "plotLR")
#' lims <- range(c(plr$xl, plr$xr, plr$yl, plr$yr), finite = TRUE)
#' g0 <- g0 + ggplot2::coord_fixed(xlim = lims, ylim = lims)
#'
#' ## Log-link fit
#' mod <- fit_log
#' theta <- mod$thetaHat
#' famly <- mod$family
#' f.eta  <- function(s) theta[1] + theta[2] * s
#' f.mu   <- function(s) famly$linkinv(f.eta(s))
#' f.zeta <- function(s) famly$linkfun(f.mu(s)) + (y - f.mu(s)) / famly$mu.eta(f.eta(s))
#'
#' g1 <- plot_ic_mapping(cbind(zl, zr), funx = f.eta, funy = f.zeta) +
#'   ggplot2::geom_abline(linetype = 2, alpha = .6) +
#'   ggplot2::labs(x = "linear predictor", y = "working response", title = NULL)
#' plr <- attr(g1, "plotLR")
#' lims <- range(c(plr$xl, plr$xr, plr$yl, plr$yr), finite = TRUE)
#' g1 <- g1 + ggplot2::coord_fixed(xlim = lims, ylim = lims)
#'
#' ggpubr::ggarrange(plotlist = list(g0, g1), nrow = 1)
#' }
#'
#' @export
plot_ic_mapping <- function (LR, funx = function(x) x, funy, lseq = 10, eps=NULL, 
                             funx_args = list(), funy_args = list(),
                             path_args = list(alpha = 0.4), point_args = list(alpha = 0.7)) 
{
  # check arguments
  if (nrow(LR) == 1L) {
    stop("`LR` contains just one interval.")
  }
  if (lseq <= 1) {
    lseq <- 2
    message("`lseq` must be at least 2; using `lseq = 2`.")
  }
  path_args <- utils::modifyList(list(alpha = 0.4), path_args)
  point_args <- utils::modifyList(list(alpha = 0.7), point_args)
  
  vals <- apply(LR, 1, function(x) seq(x[1], x[2], length.out = lseq))
  
  xs <- apply(vals, 1, function(z) do.call(funx, c(list(z), funx_args)))
  ys <- apply(vals, 1, function(z) do.call(funy, c(list(z), funy_args)))
  
  df <- data.frame(
    obs = rep(seq_len(ncol(vals)), each = nrow(vals)),
    x   = as.vector(t(xs)),
    y   = as.vector(t(ys))
  )
  
  plotLR <- data.frame(
    id = seq_len(nrow(LR)),
    xl = apply(xs, 1, min),
    xr = apply(xs, 1, max),
    yl = apply(ys, 1, min),
    yr = apply(ys, 1, max)
  )
  
  plotLR$dx <- plotLR$xr - plotLR$xl
  plotLR$dy <- plotLR$yr - plotLR$yl
  
  if (is.null(eps)) {
    
    plt <- ggplot(df, aes(x = x, y = y, group = obs)) +
      do.call(geom_path, path_args) +
      theme_bw()
    
  } else {
    
    # Draw a point when neither the horizontal nor the vertical range is large enough to be visible
    small <- plotLR$id[plotLR$dx <= eps & plotLR$dy <= eps]
    large <- setdiff(plotLR$id, small)
    
    df_lines <- df[df$obs %in% large, , drop = FALSE]
    
    df_points <- data.frame(
      obs = small,
      x = (plotLR$xl[small] + plotLR$xr[small]) / 2,
      y = (plotLR$yl[small] + plotLR$yr[small]) / 2
    )
    
    plt <- ggplot(df, aes(x = x, y = y, group = obs)) +
      do.call(geom_path, c(list(data = df_lines), path_args)) +
      do.call(geom_point, c(
        list(
          data = df_points,
          mapping = aes(x = x, y = y),
          inherit.aes = FALSE
        ),
        point_args
      )
      ) +
      theme_bw()
    
  }
  
  attr(plt, "plotLR") <- plotLR
  plt
}


utils::globalVariables(c("obs","x","y"))
