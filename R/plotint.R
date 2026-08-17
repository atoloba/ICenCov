#' Plot observed intervals
#'
#' Creates a plot from a dataframe containing left/right endpoints. 
#' Optionally, points for the true covariate \code{z} (e.g., in simulations) 
#' can be overlaid, and rows can be ordered by an outcome variable. 
#' If your column names differ, specify them through \code{zl}, \code{zr}, \code{y}, and \code{z}.
#'
#' @param dt A data frame containing the interval endpoints and, optionally,
#'   response and true covariate values.
#' @param order.by.y Logical; if \code{TRUE}, order intervals by \code{y}
#'   (ascending). If \code{FALSE}, order by the interval midpoint
#'   \eqn{(zl + zr)/2} (default).
#' @param show.z Logical; if \code{TRUE}, plot the true covariate values from
#'   the column specified by \code{z}.
#' @param zl Character string giving the name of the left-endpoint column.
#'   Defaults to \code{"zl"}.
#' @param zr Character string giving the name of the right-endpoint column.
#'   Defaults to \code{"zr"}.
#' @param y Character string giving the name of the response column, used when
#'   \code{order.by.y = TRUE}. Defaults to \code{"y"}.
#' @param z Character string giving the name of the true covariate column, used
#'   when \code{show.z = TRUE}. Defaults to \code{"z"}.
#' @param eps Optional non-negative numeric value. A point is added for intervals
#'   with width smaller than or equal to \code{eps}. If \code{NULL}, points are
#'   used only for exact observations.
#' @param interval_args Named list of additional arguments passed to
#'   \code{\link[ggplot2]{geom_linerange}}, such as \code{color},
#'   \code{linewidth}, \code{linetype}, and \code{alpha}. Defaults to
#'   \code{list(linewidth = 0.4)}.
#' @param point_args Named list of additional arguments passed to
#'   \code{\link[ggplot2]{geom_point}} for exact or sufficiently small intervals,
#'   such as \code{color}, \code{fill}, \code{size}, \code{shape},
#'   \code{alpha}, and \code{stroke}. Defaults to \code{list(size = 0.1)}.
#' @param zpoint_args Named list of additional arguments passed to
#'   \code{\link[ggplot2]{geom_point}} for the true covariate \code{z} when
#'   \code{show.z = TRUE}, such as \code{color}, \code{fill}, \code{size},
#'   \code{shape}, \code{alpha}, and \code{stroke}. Defaults to
#'   \code{list(shape = 20, color = "blue", alpha = 0.5)}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' data(carotenoids)
#' library(ggplot2)
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
plotint <- function(dt, order.by.y = FALSE, show.z = FALSE, 
                    zl = "zl", zr = "zr", y = "y", z = "z", eps = NULL,
                    interval_args = list(linewidth = 0.4),
                    point_args = list(size = 0.1),
                    zpoint_args = list(shape = 20, color = "blue", alpha = 0.5)
                    ) {
  
  if (!is.data.frame(dt))
    stop("`dt` must be dataframe")

  ncolstocheck <- c(zl, zr, if (order.by.y) y, if (show.z) z)
  if (!all(ncolstocheck %in% names(dt)))
    stop(
      "Missing required column(s): ",
      paste(setdiff(ncolstocheck, names(dt)), collapse = ", ")
    )
  
  interval_args <- utils::modifyList(list(linewidth = 0.4), interval_args)
  point_args <- utils::modifyList(list(size = 0.1), point_args)
  zpoint_args <- utils::modifyList(list(shape = 20, color = "blue", alpha = 0.5), zpoint_args)
  
  dto <- if (order.by.y) {
    dt[order(dt[[y]]), ]
  } else {
    dt[order((dt[[zl]] + dt[[zr]]) / 2), ]
  }
  
  dto$id <- 1:nrow(dto)
  
  width <- dto[[zr]] - dto[[zl]]
  
  small <- if (is.null(eps)) {
    width == 0
  } else {
    width <= eps
  }
  
  g <- ggplot(dto, aes(y = id)) +
    do.call(geom_linerange, c(
      list(mapping = aes(xmin = .data[[zl]], xmax = .data[[zr]])),
      interval_args
      )) +
    do.call(geom_point, c(
      list(data = dto[small, ], 
           mapping = aes(x = .data[[zl]])),
      point_args
      ))
  
  if (show.z) {
    g <- g +
      do.call(geom_point, c(
        list(mapping = aes(x = .data[[z]])),
        zpoint_args
      ))
    }
  
  g <- g +
    labs(y = "obs. num.", x = "Observed intervals", title = NULL) +
    theme_bw()
  
  return(g)
}

utils::globalVariables(c("id", ".data"))