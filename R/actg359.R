#' AIDS Clinical Trials Group Study 359
#'
#' The data (\eqn{n = 81}) is from a randomized clinical trial designed to 
#' compare six different antiretroviral treatment regimens in HIV-infected 
#' individuals who had previously failed the indinavir therapy, an inhibitor 
#' of the HIV (Gulick et al., 2000). One research question of interest 
#' was whether delays in initiating the new treatment could influence viral 
#' load levels at baseline, and so potentially affect treatment prognosis.
#' Indinavir failure was defined as a rebound in viral load (\eqn{> 500} copies/mL)
#' after an initial period of virologic suppression, with viral load monitored
#' at regular visits. Therefore, the time between indinavir failure and study 
#' start is interval-censored.
#' 
#' @format A data frame with 81 rows and the following variables:
#' \describe{
#'   \item{logRNA}{Log\eqn{_{10}} copies of RNA (viral load).}
#'   \item{RNA}{Copies of RNA (viral load).}
#'   \item{age}{Age (years).}
#'   \item{zl}{Time from \emph{first} viral load \eqn{> 500} copies/mL to randomization, in weeks.}
#'   \item{zr}{Time from \emph{last} viral load \eqn{< 500} copies/mL to randomization, in weeks.}
#' }
#'
#' @details
#'
#' Data may be used in publications that cite the references below. 
#'
#'
#' @references
#' Gulick, R. M., Hu, X. J., Fiscus, S. A., Fletcher, C. V., Haubrich, R., Cheng, H., Acosta, E.,
#' Lagakos, S. W., Swanstrom, R., Freimuth, W., Snyder, S., Mills, C., Fischl, M., Pettinelli, C.,
#' & Katzenstein, D. (2000). Randomized study of saquinavir with ritonavir or nelfinavir together
#' with delavirdine, adefovir, or both in human immunodeficiency virus-infected adults with
#' virologic failure on indinavir: AIDS Clinical Trials Group Study 359.
#' \emph{The Journal of Infectious Diseases}, 182(5), 1375–1384. \doi{10.1086/315867}
#'
#' @usage data(actg359)
#' @docType data
#' @keywords datasets
#' @name actg359
#' @title AIDS Clinical Trials Group Study 359
NULL