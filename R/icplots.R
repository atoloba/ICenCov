# Has geom_rect option


#' Plot observed intervals
#'
#' Creates a plot from a dataframe containing left/right endpoints \code{zl},
#' \code{zr}. Optionally, points for the true covariate \code{z}
#' (e.g., in simulations) can be overlaid, and rows can be ordered by
#' an outcome \code{y}.
#'
#' @param dt A data frame with columns \code{zl}, \code{zr}. Optional columns:
#'   \code{y} (used only when \code{order.by.y = TRUE}) and \code{z}
#'   (used only when \code{show.z = TRUE}).
#' @param show.z Logical; if \code{TRUE}, plot the true \code{z} in blue.
#' @param order.by.y Logical; if \code{TRUE}, order intervals by \code{y}
#'   (ascending). If \code{FALSE}, order by the interval midpoint
#'   \eqn{(zl + zr)/2} (default).
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' data(carotenoids)
#'
#' plotint(carotenoids) +
#'  labs(x=expression(paste("Total carotenoid concentration [", mu, "mol/L]")))
#'
#'
#' data(actg359)
#'
#' plotint(actg359) +
#'  labs(x="Time since indinavir failure (weeks)")
#'
#'
#' @importFrom ggplot2 ggplot aes geom_linerange geom_point labs theme_bw
#' @export
plotint <- function(dt, order.by.y = FALSE, show.z = FALSE){
  if(!is.data.frame(dt)) stop("`dt` must be dataframe")
  if(order.by.y){
    dto <- dt[order(dt$y),]
  } else {
    dto <- dt[order((dt$zl+dt$zr)/2),]
  }
  dto$id <- 1:nrow(dto)

  g <- ggplot(dto, aes(y=id)) +
    geom_linerange(aes(xmin = zl, xmax = zr), linewidth=.4)
  if(show.z)
    g <- g + geom_point(aes(x=z), shape=20, color="blue", alpha=.5)
  g <- g + labs(y = "obs. num.", x = "Observed intervals", title="") +
    theme_bw()
  return(g)
}




#' Plot augmented Turnbull estimator(s)
#'
#' Plots the augmented Turnbull estimator returned by \code{icglm}. Accepts a list
#' \code{lwHat} of one or more estimators for the same censored covariate and overlays them.
#'
#' @param lwHat List of augmented Turnbull estimators returned by \code{icglm}.
#' @param col Vector of colors (one per element in \code{lwHat}).
#' @param nleg Character string used as legend title.
#' @param gini Optional \pkg{ggplot2} object to add layers to; if \code{NULL},
#' a new plot is created.
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' data(actg359)
#' linmod <- icglm(logRNA ~ age + ic(zl,zr,"waitime"), family = gaussian,
#'                 data = actg359, Lin = F, Rin = T)
#' gammod <- icglm(RNA ~ age + ic(zl,zr,"waitime"), family = Gamma("log"),
#'                 data = actg359, Lin = F, Rin = T)
#' plot.wHat(list("LM" = linmod$wHat, "GLM" = gammod$wHat))
#'
#' @importFrom ggplot2 ggplot aes geom_linerange geom_segment labs scale_color_manual theme_bw
#' @export
plot.wHat <- function(lwHat, col = 1:length(lwHat), nleg = "",
                      gini = NULL){
  if(length(col) != length(lwHat)) stop("Provide one color per wHat element")
  if(is.null(names(lwHat))) names(lwHat) <- "wHat"
  names(col) <- names(lwHat)

  if(!is.null(gini)) g <- gini else g <- ggplot()

  for(i in 1:length(lwHat)){
    omega <- lwHat[[i]][,-1]
    class(omega) <- "numeric"
    preplot1a <- data.frame(a = c(0, omega[,2]),
                            b = c(omega[,1], omega[nrow(omega),2]),
                            pfun = c(0, omega[,4]),
                            group = names(lwHat)[i])
    preplot1b <- data.frame(x1 = omega[,1],
                            x2 = omega[,2],
                            y1 = c(0, omega[-nrow(omega), 4]),
                            y2 = omega[,4],
                            group = names(lwHat)[i])
    g <- g +
      geom_linerange(aes(xmin = a, xmax = b, y = pfun, color=group), alpha=.5, preplot1a) +
      geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), colour=col[i], alpha=.5, preplot1b)
      # geom_rect(aes(xmin = x1, xmax = x2, ymax = y1, ymin = y2), fill=col[i], colour=alpha(col[i],.3), alpha=.5, preplot1b)
  }
  g +
    labs(x = "z", y = expression(hat(W))) +
    scale_color_manual(name = nleg, values = col) +
    theme_bw() +
    theme(legend.position = "right")
}
