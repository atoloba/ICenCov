#' Study of carotenoids
#'
#' Baseline visit data for a subgroup (\eqn{n = 230}) of participants in the
#' PREDIMED-Plus study, a Spanish primary-prevention trial of cardiovascular
#' disease in older adults (55–75 years) at high cardiovascular risk.
#'
#' @format A data frame with 230 rows and the following variables:
#' \describe{
#'   \item{id}{Participant identifier.}
#'   \item{zl, zr}{Left and right endpoints of circulating carotenoids (interval-censored).}
#'   \item{acarl, acarr}{Left and right endpoints of \eqn{\alpha}-carotene (interval-censored).}
#'   \item{age}{Age in years.}
#'   \item{sex}{Factor with 2 levels: \code{"Male"}, \code{"Female"}.}
#'   \item{bmi}{Body mass index (kg/m\eqn{^2}).}
#'   \item{energy}{Total energy intake (Mcal/day).}
#'   \item{glucose}{Fasting glucose level (mg/dL).}
#'   \item{diabetes}{Indicator of diabetes: treated diabetes, HbA1c \eqn{\ge} 6.5,
#'     or fasting glucose \eqn{\ge} 126 mg/dL.}
#' }
#'
#' @details
#' See the PREDIMED-Plus website for study information and protocols.
#' Data may be used in publications that cite the reference below.
#'
#' @source PREDIMED-Plus Study. \url{https://www.predimedplus.com/en/}
#'
#' @references
#' Marhuenda-Muñoz, M., Domínguez-López, I., Langohr, K., Tresserra-Rimbau, A., ... 
#' Gómez Melis, G., Lamuela-Raventós, R. M. (2022).
#' Circulating carotenoids are associated with favorable lipid and fatty acid profiles in an older
#' population at high cardiovascular risk. \emph{Frontiers in Nutrition}, 9, 967967.
#' \doi{10.3389/fnut.2022.967967}
#'
#' @usage data(carotenoids)
#' @docType data
#' @keywords datasets
#' @name carotenoids
#' @title Study of carotenoids
NULL





