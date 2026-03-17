#' Pearson, deviance and quantile residuals for a GLM with an interval-censored covariate
#'
#' Computes residuals for a GLM
#' \eqn{\mu=\mathrm{linkinv}(X^\top\beta + Z\gamma)},
#' where \code{Z} is only known to lie between \code{zl} and \code{zr}. Residuals are returned
#' in two ways: (i) as functions \eqn{r(z)} for a supplied value of \code{z} (via \code{fpear},
#' \code{fdev}, \code{fquant}), and (ii) as interval-censored residuals with endpoints
#' \code{rl} and \code{rr} (via \code{icpear}, \code{icdev}, \code{icquant}).
#'
#'
#' @param theta vector of parameter estimates \eqn{(\hat\beta,\hat\gamma,\hat\phi)};
#'    \eqn{\hat\phi} is only used for quantile residuals.
#' @param y Response vector.
#' @param zl,zr Numeric vectors with left and right endpoints of the
#'   observed intervals for the interval-censored covariate.
#' @param X Design matrix for the non-censored covariates.
#' @param famly A family object returned by \code{\link{icglm}}.
#' @param type Character vector of residual types (\code{"pearson"}, \code{"deviance"},
#'   \code{"quantile"}); \code{"all"} (default).
#'
#' @return A list.
#'
#' @details
#' Quantile residuals are based on the Probability Integral Transform (PIT).
#' For continuous families the PIT \eqn{u=F(y\mid\hat\mu)} is uniquely determined by \eqn{\hat\mu},
#' and the quantile residual is \eqn{r=\Phi^{-1}(u)}, where \eqn{\Phi} is the standard normal cdf.
#' When \eqn{\hat\mu} is interval-censored with endpoints \code{mul} and \code{mur}, the PIT
#' is between \eqn{F(y\mid \mathrm{mur})} and \eqn{F(y\mid \mathrm{mul})}.
#' Function \code{fquant} computes \eqn{r=\Phi^{-1}(u)} at a supplied \eqn{\hat\mu},
#' while \code{icquant} corresponds to the residual interval induced by \code{mul} and \code{mur}.
#'
#' For discrete families (\code{binomial} and \code{poisson}), given \eqn{\hat\mu} the PIT is an interval
#' \eqn{u\in[F(y^-\mid\hat\mu),\,F(y\mid\hat\mu)]} due to jumps in the cdf. With interval-censored \eqn{\hat\mu},
#' these endpoints are combined with \code{mul} and \code{mur} to form \code{icquant}.
#' In this case, \code{fquant} returns *randomized* quantile residuals by sampling
#' \eqn{u^* \sim \mathrm{Unif}(F(y^-\mid\hat\mu),\,F(y\mid\hat\mu))} and returning \eqn{\Phi^{-1}(u^*)}.
#'
#'
#' @seealso \code{\link{icglm}}
#'
#' @examples
#' data(actg359)
#'
#' fit2 <- icglm(RNA ~ age + ic(zl, zr, "waitime"), family = Gamma("log"),
#'             data = actg359, Lin = TRUE, Rin = FALSE)
#' out <- ires(fit2$thetaHat, actg359$RNA, actg359$zl, actg359$zr,
#'             model.matrix(RNA~age, actg359), fit2$family)
#'
#' @importFrom stats qnorm runif
#' @export
ires <- function(theta, y, zl, zr, X, famly, type="all") {

  if("all" %in% type) type <- c("pearson","deviance","quantile")
  if(!all(type %in% c("pearson","deviance","quantile")))
    stop("type should be one of 'deviance', 'pearson', 'quantile'")

  y  <- as.numeric(y)
  zl <- as.numeric(zl)
  zr <- as.numeric(zr)

  mu <- funmu(theta, zl, zr, X, famly$linkinv)
  if(!is.null(mu$error)) return(mu)
  out <- list(fmu=mu$fmu, icmu=mu$icmu)
  fmu <- mu$fmu
  mul <- mu$icmu$mul
  mur <- mu$icmu$mur


  if("pearson"%in%type)
    # typically decreasing except for Tweedie with p>2,
    # for which it's u-shaped with the min at p*y/(p-2)
  {
    fpear <- function(val, ismu = FALSE) {
      if (!(length(val) %in% c(1, nrow(X))))
        stop(paste("input length must be 1 or", nrow(X)))

      mu <- if (ismu) val else fmu(val)

      if (famly$family == "inverse.gaussian" && any(mu <= 0, na.rm = TRUE))
        stop("inverse.gaussian: mu must be > 0")

      (y - mu) / sqrt(famly$variance(mu))
    }

    rl <- (y - mur) / sqrt(famly$variance(mur))
    rr <- (y - mul) / sqrt(famly$variance(mul))

    # Inverse gaussian (Tweedie with p=3)
    if(famly$family=="inverse.gaussian"){
      muleps <- pmax(mul, .Machine$double.eps^0.5)
      mureps <- pmax(mur, .Machine$double.eps^0.5)
      rl0 <- (y - mureps) / sqrt(famly$variance(mureps))
      rr0 <- (y - muleps) / sqrt(famly$variance(muleps))
      rr <- pmax(rl0, rr0)
      rl <- ifelse(muleps <= 3*y & 3*y <= mureps, fpear(3*y, ismu = TRUE), pmin(rl0, rr0))
    }

    icpear <- data.frame(rl,rr)
    out$fpear  <- fpear
    out$icpear <- icpear
  }


  if("deviance"%in%type)
    # Monotone decreasing in mu, equals zero at mu=y
  {
    wt <- rep(1, length(y))

    fdev <- function(val, ismu = FALSE) {
      if (!(length(val) %in% c(1, nrow(X))))
        stop(paste("input length must be 1 or", nrow(X)))
      mu <- if (ismu) val else fmu(val)
      sign(y - mu) * sqrt(famly$dev.resids(y, mu, wt = wt))
    }

    rl <- sign(y - mur) * sqrt(famly$dev.resids(y, mur, wt = wt))
    rr <- sign(y - mul) * sqrt(famly$dev.resids(y, mul, wt = wt))

    icdev <- data.frame(rl = rl, rr = rr)
    out$fdev  <- fdev
    out$icdev <- icdev
  }


  if ("quantile" %in% type)  {
    phi  <- theta[length(theta)]
    eps <- .Machine$double.eps^0.5

    if (famly$family %in% c("poisson", "binomial")) {
      # Discrete case:
      # For fixed mu, U = F(y | mu) interval is [F(y^-), F(y)] = [F(y) - P(Y=y), F(y)].
      # Over mu in [mul, mur] and with F decreasing in mu:
      #   lower U endpoint = F(y^- | mu = mur)   (smallest, since mu largest)
      #   upper U endpoint = F(y   | mu = mul)   (largest, since mu smallest)

      F_mur <- pfam(y, phi = phi, famly = famly, mu = mur)
      P_mur <- dfam(y, phi = phi, famly = famly, mu = mur)
      F0_mur <- pmax(F_mur - P_mur, 0)  # left-limit F(y^-)

      F_mul <- pfam(y, phi = phi, famly = famly, mu = mul)

      # clamp away from 0/1 to avoid +/-Inf
      uL <- pmin(pmax(F0_mur, eps), 1 - eps)
      uR <- pmin(pmax(F_mul,  eps), 1 - eps)

      rl <- qnorm(uL)
      rr <- qnorm(uR)

    } else {
      # Continuous case:
      # U = F(y | mu). Over mu in [mul, mur] and F decreasing in mu:
      # U ranges from F(y | mur) (smallest) to F(y | mul) (largest)

      F_mur <- pfam(y, phi = phi, famly = famly, mu = mur)
      F_mul <- pfam(y, phi = phi, famly = famly, mu = mul)

      uL <- pmin(pmax(F_mur, eps), 1 - eps)
      uR <- pmin(pmax(F_mul, eps), 1 - eps)

      rl <- qnorm(uL)
      rr <- qnorm(uR)
    }

    fquant <- function(val, ismu = FALSE) {
      # For discrete families returns *randomized* quantile residuals
      # qnorm(U), U ~ Unif(F(y^-), F(y))

      if (!(length(val) %in% c(1, nrow(X))))
        stop(paste("input length must be 1 or", nrow(X)))

      mu <- if (ismu) val else fmu(val)

      # cdf at observed y
      Fy <- pfam(y, phi = phi, famly = famly, mu = mu)

      if (famly$family %in% c("poisson", "binomial")) {
        # pmf at observed y
        Py <- dfam(y, phi = phi, famly = famly, mu = mu)

        # F(y^-) = F(y) - P(Y=y)
        F0 <- pmax(Fy - Py, 0)

        # clamp away from 0/1 to avoid +/-Inf
        Fy  <- pmin(pmax(Fy,  eps), 1 - eps)
        F0 <- pmin(pmax(F0, eps), 1 - eps)

        # Randomized PIT
        u <- runif(length(Fy), min = F0, max = Fy)

        return(qnorm(u))
      } else {
        # Continuous: unique PIT
        u <- pmin(pmax(Fy, eps), 1 - eps)
        return(qnorm(u))
      }
    }

    out$fquant <- fquant
    out$icquant <- data.frame(rl,rr)
  }

  return(out)
}



#' Fitted mean for a GLM with an interval-censored covariate
#'
#' Computes the fitted mean of a GLM
#' \eqn{\mu=\mathrm{linkinv}(X^\top\beta + Z\gamma)},
#' where \code{Z} is only known to lie between
#' \code{zl} and \code{zr}, in two ways: (i) as a function \eqn{\hat\mu(z)} for a supplied value of
#' \code{z}, and (ii) as an interval-censored fitted mean with endpoints \code{mul}
#' and \code{mur}.
#'
#' @param theta Numeric vector of parameter estimates \eqn{(\hat\beta,\hat\gamma,\hat\phi)};
#'   \eqn{\hat\phi} is not used here.
#' @param zl,zr Numeric vectors with left and right endpoints of the
#'   observed intervals for the interval-censored covariate.
#' @param X Design matrix for the non-censored covariates.
#' @param linkinv Inverse link function of the GLM.
#'
#' @return A list with components:
#' \describe{
#'   \item{fmu}{A function that evaluates the fitted mean for a supplied value of the interval-censored covariate.}
#'   \item{icmu}{A data frame with columns \code{mul} and \code{mur} giving the endpoints of the fitted mean for each observation.}
#'   \item{error}{\code{NULL} if no inconsistencies are detected.}
#' }
#'
#' @export
funmu <- function(theta, zl, zr, X, linkinv){

  zl <- as.numeric(zl)
  zr <- as.numeric(zr)

  lbet <- length(theta) - 2
  mul <- linkinv(X%*%theta[1:lbet] + theta[lbet+1]*zl)
  mur <- linkinv(X%*%theta[1:lbet] + theta[lbet+1]*zr)

  mudif <- mur - mul
  tol <- .Machine$double.eps^0.5
  allok   <- all(mudif >= -tol)
  allrev  <- all(mudif <=  tol)

  if (allrev && !allok){
    tmp <- mul; mul <- mur; mur <- tmp
  }

  return(list(
    fmu = function(val){
      if(!(length(val) %in% c(1, nrow(X))))
        stop(paste("input length must be 1 or", nrow(X)))
      return(linkinv(X%*%theta[1:lbet] + theta[lbet+1]*val))
    },
    icmu = data.frame(mul, mur),
    error = if(!allrev && !allok) "Inconsistent monotonicity: for some rows mu(zr) < mu(zl) and for others mu(zr) > mu(zl)" else NULL
  ))
}










