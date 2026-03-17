#' Plot observed intervals
#'
#' Creates a plot from a dataframe containing left/right endpoints (by default
#' \code{zl} and \code{zr}). Optionally, points for the true covariate \code{z}
#' (e.g., in simulations) can be overlaid, and rows can be ordered by
#' an outcome \code{y}. If your column names differ, pass them in \code{namecols}.
#'
#' @param dt A data frame with columns \code{zl}, \code{zr}. Optional columns:
#'   \code{y} (used only when \code{order.by.y = TRUE}) and \code{z}
#'   (used only when \code{show.z = TRUE}).
#' @param show.z Logical; if \code{TRUE}, plot the true \code{z} in blue.
#' @param order.by.y Logical; if \code{TRUE}, order intervals by \code{y}
#'   (ascending). If \code{FALSE}, order by the interval midpoint
#'   \eqn{(zl + zr)/2} (default).
#' @param namecols Named list giving the column names in \code{dt} that correspond to
#'   the interval endpoints \code{zl} and \code{zr}, and (if used) \code{z} and \code{y}.
#' @param eps Non-negative numeric tolerance used to detect intervals that should
#'   be treated as points for plotting, i.e. when \code{zl} and \code{zr} are
#'   equal up to numerical precision.
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
#' @importFrom utils modifyList
#' @export
plotint <- function (dt, order.by.y = FALSE, show.z = FALSE,
                     namecols = list(zl = "zl", zr = "zr", y = "y", z = "z"),
                     eps = .Machine$double.eps^0.5) {

  if (!is.data.frame(dt))
    stop("`dt` must be dataframe")

  # fill missing names with defaults
  ncols <- namecols
  ncols <- modifyList(list(zl = "zl", zr = "zr", y = "y", z = "z"), ncols)
  ncolstocheck <- c(ncols$zl, ncols$zr, if (order.by.y) ncols$y, if (show.z) ncols$z)
  if (!all(ncolstocheck %in% names(dt)))
    stop("Missing required column(s): ", paste(setdiff(ncolstocheck, names(dt)), collapse = ", "))

  dto <- if (order.by.y) dt[order(dt[[ncols$y]]), ]
         else dt[order((dt[[ncols$zl]] + dt[[ncols$zr]])/2), ]
  dto$id <- 1:nrow(dto)

  g <- ggplot(dto, aes(y = id)) +
    geom_linerange(aes(xmin = .data[[ncols$zl]], xmax = .data[[ncols$zr]]), linewidth = 0.4) +
    geom_point(data = dto[(dto[[ncols$zl]] - dto[[ncols$zr]]) > -eps, ],
               aes(x = .data[[ncols$zl]]), size= 0.1)

  if (show.z)
    g <- g + geom_point(aes(x = .data[[ncols$z]]),
                        shape = 20, color = "blue", alpha = 0.5)
  g <- g + labs(y = "obs. num.", x = "Observed intervals",
                title = NULL) + theme_bw()
  return(g)
}

utils::globalVariables(c("id", ".data"))
