#' Plot (augmented) Turnbull estimator(s)
#'
#' Plots the (augmented) Turnbull cumulative distribution function (CDF) for one
#' interval-censored covariate. When several models are fit with the same
#' interval-censored covariate (e.g., different model covariates),
#' the corresponding CDFs can be overlaid. The CDF can be drawn
#' as piecewise segments (augmented Turnbull) or as Turnbull-style blocks.
#'
#' @param lwhat List of (augmented) Turnbull estimates in the format returned by \code{\link{icglm}}.
#' @param col Vector of colors (one per element in \code{lwhat}).
#' @param nleg Character string used as legend title.
#' @param style Character vector. Either a single value (applied to all) or one value per element of \code{lwhat}:
#'   \code{"seg"} for piecewise segments (augmented Turnbull) and \code{"block"} for Turnbull-style blocks.
#' @param gini Optional \pkg{ggplot2} object to add layers to; if \code{NULL},
#'   a new plot is created.
#' @param line_args Named list of arguments passed to
#'   \code{\link[ggplot2]{geom_linerange}}, such as \code{linewidth},
#'   \code{linetype}, and \code{alpha}. Color is controlled by \code{col}
#'   and should not be supplied. Defaults to \code{list(alpha = 0.5)}.
#' @param segment_args Named list of arguments for \code{style = "seg"}, passed to
#'   \code{\link[ggplot2]{geom_segment}}, such as \code{linewidth},
#'   \code{linetype}, and \code{alpha}. Color is controlled by \code{col}
#'   and should not be supplied. Defaults to \code{list(alpha = 0.5)}.
#' @param block_args Named list of arguments for \code{style = "block"}, passed to
#'   \code{\link[ggplot2]{geom_rect}}, such as \code{linewidth},
#'   \code{linetype}, and \code{alpha}. Color and fill are controlled by
#'   \code{col} and should not be supplied. Defaults to
#'   \code{list(alpha = 0.5)}.
#'   
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' data(actg359)
#' linmod <- icglm(logRNA ~ age + ic(zl,zr,"waitime"), family = gaussian,
#'                 data = actg359, Lin = TRUE, Rin = FALSE)
#' gammod <- icglm(RNA ~ age + ic(zl,zr,"waitime"), family = Gamma("log"),
#'                 data = actg359, Lin = TRUE, Rin = FALSE)
#' plot_ic_cdf(list("LM" = linmod$wHat, "GLM" = gammod$wHat))
#'
#' @importFrom ggplot2 ggplot aes geom_linerange geom_segment geom_rect geom_blank
#' @importFrom ggplot2 labs scale_color_manual theme_bw theme
#' @export
plot_ic_cdf <- function(lwhat, col = seq_along(lwhat), nleg = "",
                        style = "seg", gini = NULL,
                        line_args = list(alpha = 0.5),
                        segment_args = list(alpha = 0.5),
                        block_args = list(alpha = 0.5)
                        ) {

  if(length(col) != length(lwhat)) stop("Provide one color per wHat element")
  if (is.null(names(lwhat))) {
    names(lwhat) <- if (length(lwhat) == 1L) {
      "wHat"
    } else {
      paste0("wHat", seq_along(lwhat))
    }
  }
  names(col) <- names(lwhat)

  if (length(style) == 1) style <- rep(style, length(lwhat))
  if (length(style) != length(lwhat)) stop("`style` must have length 1 or length(lwhat)")
  style <- match.arg(style, choices = c("seg", "block"), several.ok = TRUE)

  line_args <- utils::modifyList(list(alpha = 0.5), line_args)
  segment_args <- utils::modifyList(list(alpha = 0.5), segment_args)
  block_args <- utils::modifyList(list(alpha = 0.5), block_args)
  
  g <- if (!is.null(gini)) gini else ggplot()


  # Preset the correct x-limits for later layers (e.g., stat_function in plot_icres_cdf function)
  x_rng <- range(unlist(lapply(lwhat, function(w) {
    omega <- w[, -1, drop = FALSE]
    c(as.numeric(omega[, 1]), as.numeric(omega[, 2]))
  })), finite = TRUE)
  g <- g + geom_blank(data = data.frame(x = x_rng), aes(x = x, y = 0))

  # Small numerical padding at the limits
  pad <- .Machine$double.eps^0.5
  
  # Add the ic_cdf layers
  for(i in seq_along(lwhat)){
    
    omega <- lwhat[[i]][, -1, drop = FALSE]
    storage.mode(omega) <- "numeric"
    preplotA <- data.frame(a = c(min(omega[,1])-pad, omega[,2]),
                           b = c(omega[,1], max(omega[,2])+pad),
                           pfun = c(0, omega[,4]),
                           group = names(lwhat)[i])
    preplotB <- data.frame(x1 = omega[,1],
                           x2 = omega[,2],
                           y1 = c(0, omega[-nrow(omega), 4]),
                           y2 = omega[,4],
                           group = names(lwhat)[i])
    g <- g +
      do.call(geom_linerange, c(
        list(mapping = aes(xmin = a, xmax = b, y = pfun, color = group),
             data = preplotA),
        line_args
        ))

    if (style[i] == "seg") {
      g <- g + 
        do.call(geom_segment, c(
          list(mapping = aes(x = x1, xend = x2, y = y1, yend = y2),
               color = col[i],
               data = preplotB),
          segment_args
          ))
    } else {
      g <- g + 
        do.call(geom_rect, c(
          list(mapping = aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2),
               fill  = col[i], color = scales::alpha(col[i], 0.3),
               data = preplotB),
          block_args
        ))
    }
  }

  g +
    labs(x = "z", y = expression(hat(W)), title = NULL) +
    scale_color_manual(name = nleg, values = col) +
    theme_bw() +
    theme(legend.position = "right")
}

utils::globalVariables(c("a","b","pfun","group","x","x1","x2","y1","y2"))
