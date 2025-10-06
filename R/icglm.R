#' Generalized Linear Models with an Interval-Censored Covariate
#'
#' Fits a GLM when one covariate is only known to fall within an interval.
#' It uses the **GELc** estimator: a semiparametric, likelihood-based approach
#' built on an augmented version of Turnbull's nonparametric estimator for
#' interval-censored data.
#'
#' @param formula A model formula of the form
#'   `y ~ x1 + x2 + ... + ic(zl, zr, "z")`.
#'   The interval-censored covariate is supplied via \code{ic()}, where
#'   \code{zl} and \code{zr} are the left and right endpoints, and the third
#'   argument gives the covariate name used in the linear predictor.
#'   Use \code{Lin} and \code{Rin} to state whether endpoints are included
#'   in the interval (see below).
#' @param family A \code{stats::family} object. Supported families are
#'   \code{gaussian()}, \code{Gamma()}, \code{inverse.gaussian()},
#'   \code{binomial()}, and \code{poisson()}.
#' @param Lin,Rin Logical scalar or length-\eqn{n} vector indicating whether the
#'   left/right endpoints are included in the interval. Defaults are \code{TRUE} for both.
#' @param tolpar Relative change threshold for the parameter vector used as a
#'   stopping rule (default \code{1e-4}).
#' @param tolik Relative change threshold for the log-likelihood used as a
#'   stopping rule (default \code{1e-4}).
#' @param maxiter Maximum number of GELc iterations (default \code{50}).
#' @param maxsec Maximum running time in seconds. By default,
#' \code{maxsec = maxh * 60^2}.
#' @param maxh Maximum running time in hours (default \code{12}).
#' @param start Optional numeric vector of length \code{p + 2} giving initial
#'   values for \code{c(beta, gamma, phi)}, where \code{p} is the length of beta.
#'   If \code{NULL}, starting values are taken from
#'   \code{glm(y ~ x1 + x2 + ... + midz)} with \code{midz = 0.5 * (zl + zr)}.
#' @param PRINT Integer verbosity flag; \code{0} is silent (default),
#'   \code{1} prints iteration summaries, and \code{2} also compares analytic
#'   and numeric gradients for stability checks.
#' @param size Optional parameter for the \code{binomial} family. If \code{NULL}
#'   (default), Bernoulli trials are assumed (\code{size = 1}).
#'
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{thetaHat}}{Numeric vector of estimated parameters
#'     \eqn{(\hat\beta, \hat\gamma, \hat\phi)}.}
#'   \item{\code{wHat}}{Matrix with the estimated density over the
#'     augmented Turnbull intervals.}
#'   \item{\code{A}}{An \eqn{n \times m} indicator matrix mapping each observed
#'     interval to the \eqn{m} augmented Turnbull intervals.}
#'   \item{\code{vcov}}{Approximate variance–covariance matrix based on the
#'     observed Hessian at \code{(thetaHat, wHat)}}
#'   \item{\code{convergence}}{A list with:
#'   \itemize{
#'          \item{\code{convOPT, mssOPT}}: integer code and message returned by the optimizer.
#'          \item \code{reachmaxi, reachmaxt}: flags for \code{maxiter} or \code{maxtime} were reached.
#'          \item \code{reachtolik}: \code{TRUE} if convergence was achieved in a plateau of the loglikelihood function.
#'          \item \code{iter_vals, iter_reldif}: per-iteration parameter values and relative diffs.
#'          \item \code{start_used}: starting values used.
#'    }
#'   \item{\code{call,family}}{Info.}
#' }
#'
#'
#' @examples
#' ## --- Logistic regression -----------------------------------------
#' data(carotenoids)
#' caro <- subset(carotenoids, !is.na(diabetes))
#'
#' fit1 <- icglm(diabetes ~ sex + bmi + ic(zl, zr, "caro"), family = binomial,
#'   data = caro, Lin = TRUE,
#'   Rin  = with(caro, zl == zr)
#' )
#' summary(fit1)
#'
#' # Odds ratio (OR) and 95% CI for 'caro':
#' i = 4
#' OR_caro  <- exp(coef(fit1)[i])
#' CI_caro  <- exp(coef(fit1)[i] + c(-1,1)*qnorm(0.975)*sqrt(vcov(fit1)[i,i]))
#' cat(sprintf("OR = %.2f (%.2f, %.2f)\n",
#'     OR_caro, CI_caro[1], CI_caro[2]))
#'
#'
#' ## --- Continuous outcome (Gamma with log link) --------------------------
#' data(actg359)
#'
#' fit2 <- icglm(RNA ~ age + ic(zl, zr, "waitime"), family = Gamma("log"),
#'   data = actg359, Lin = FALSE, Rin = TRUE)
#' summary(fit2)
#'
#' # Multiplicative effect of a 3-week increase in 'waitime' on the response scale:
#' i = 3
#' mult_wt   <- exp(3 * coef(fit2)[i]) - 1
#' CI_wt <- exp(3 * coef(fit2)[i] + 3*c(-1,1)*qnorm(0.975)*sqrt(vcov(fit2)[i,i])) - 1
#' cat(sprintf("OR = %.1f (%.1f, %.1f)\n",
#'     mult_wt*100, CI_wt[1]*100, CI_wt[2]*100))
#'
#'
#' @seealso \code{\link[stats]{glm}},
#'
#' @references
#' Gómez, G., Espinal, A., & Lagakos, S. W. (2003).
#' Inference for a linear regression model with an interval-censored covariate.
#' \emph{Statistics in Medicine}, 22(3), 409–425. \doi{10.1002/sim.1326}
#'
#' @export
icglm <- function(formula, family = gaussian, data,
                  Lin = TRUE, Rin = TRUE,
                  tolpar=1e-4, tolik=1e-4,
                  maxiter = 50, maxsec = maxh*60^2, maxh = 12,
                  start = NULL, PRINT = 0, size = NULL){

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }

  if (missing(data)) data <- environment(formula)

  # Extract elements y, zl, zr, X
  tf <- terms(formula, specials = "ic", data = data)
  sp <- attr(tf, "specials")$ic
  if(length(sp) < 1)
    stop("Not detected ic term")
  if (length(sp) > 1)
    stop("GELc currently not available for >1 ic term")

  mf <- model.frame(tf, data = data)
  y <- unclass(model.response(mf))
  X <- model.matrix(drop.terms(tf, sp-1, keep.response = FALSE), mf)
  zl <- mf[[sp]][,1]
  zr <- mf[[sp]][,2]

  spt <- colnames(mf[[sp]])[1]
  attr(zl,"name") <- substr(spt, 1, nchar(spt)-1)

  # GELc estimator
  fwd_names <- setdiff(names(formals()), c("formula", "family", "data"))
  fwd_vals  <- mget(fwd_names, inherits = TRUE)
  args <- c(list(y = y, zl = zl, zr = zr, X = X, famly = family), fwd_vals)
  fit  <- do.call(gelc.fit, args, envir = environment())

  # Output
  out <- fit[c("thetaHat","wHat","A","vcov")]
  out$convergence <- fit[c("convOPT","mssOPT","reachmaxi","reachmaxt","reachtolik",
                           "iter_vals","iter_reldif","start_used")]
  out$call <- fit$call
  out$family <- family

  class(out) <- "icglm"
  return(out)
}



gelc.fit <- function(y, zl, zr, X, famly = Gamma("log"),
                     Lin = TRUE, Rin = TRUE,
                     gr_linkinv = NULL, gr_mu.eta = NULL,
                     tolpar=1e-4, tolik=1e-4,
                     maxiter = 50, maxsec = maxh*60^2, maxh = 12,
                     start = NULL, PRINT = 0, size = NULL){

  if(is.null(gr_linkinv)|is.null(gr_mu.eta)){
  if(famly$link == "log"){
    gr_linkinv = function(eta) exp(eta)
    gr_mu.eta = function(eta) exp(eta)
  } else if(famly$link == "logit"){
    gr_linkinv = function(eta) 1 / (1 + exp(-eta))
    gr_mu.eta = function(eta) gr_linkinv(eta)*(1-gr_linkinv(eta))
  } else if(famly$link == "probit"){
    gr_linkinv = function(eta) pnorm(eta)
    gr_mu.eta = function(eta) dnorm(eta)
  } else if(famly$link == "cloglog"){
    gr_linkinv = function(eta) -expm1(-exp(eta))
    gr_mu.eta = function(eta){
      exp(eta - exp(eta))
    }
  } else if(famly$link == "cauchit"){
    gr_linkinv = function(eta) pcauchy(eta)
    gr_mu.eta = function(eta) dcauchy(eta)
  } else {
    gr_linkinv = famly$linkinv
    gr_mu.eta = famly$mu.eta
    if(!famly$link%in%c("inverse","1/mu^2","sqrt","identity"))
      warning("Default linkinv and mu.eta are taken for the analytic gradient, please check numerical stability with PRINT=2.")
  }
  }


  # Kappa matrix
  Aintmap <- Aintmap2(zl,zr, Lin, Rin)
  Kappas <- Aintmap$A
  intmap <- Aintmap$intmap
  n <- nrow(Kappas)
  m <- ncol(Kappas)
  Ipoint <- abs(intmap[2,] - intmap[1,]) < .Machine$double.eps^0.5

  zname <- if(!is.null(attr(zl,"name"))) attr(zl,"name") else "z"

  # Initial conditions
  wHat <- (Ipoint + (1-Ipoint)*(intmap[2,]-intmap[1,])) / sum(Ipoint + (1-Ipoint)*(intmap[2,]-intmap[1,]))
  if(is.null(start)){
    midz <- 0.5*(zl+zr)
    glm0 <- try(glm(y~., data = data.frame(y, X[,-1], midz), famly), silent = TRUE)
    names(glm0$coef)[1:ncol(X)] <- colnames(X)
    names(glm0$coef)[length(glm0$coef)] <- zname
    betaHat <- glm0$coef[1:(length(glm0$coef)-1)]
    gammaHat <- glm0$coef[length(glm0$coef)]
    phiHat <- summary(glm0)$dispersion
    names(phiHat) <- "disp"
  } else {
    betaHat = start[1:(length(start)-2)]
    gammaHat = start[length(start)-1]
    phiHat = start[length(start)]
    lbet <- length(betaHat)
    if(lbet != ncol(X)) stop(paste("Start is not vector of length",ncol(X)+2))
    names(betaHat) <- colnames(X)
    names(gammaHat) <- zname
    names(phiHat) <- "disp"
  }
  thetHat = c(betaHat, gammaHat, phiHat)
  lbet <- length(betaHat)
  start_used <- thetHat

  # Initialize
  lreldif_par <- vector()
  lreldif_lik <- vector()
  lthetHat <- vector()
  llikv <- vector()
  cnv <- vector()
  cnvmss <- vector()
  likv = 0
  iter = 0
  imaxi = 0
  itolik = 0
  imaxt = 0
  gotime = Sys.time()

  repeat{

    # save previous
    wOld.out <- wHat
    thetOld.out <- thetHat
    likvOld <- likv
    iter <- iter+1
    if(PRINT!=0) cat("Iter.",iter,"\n=======", fill=T)

    ## Update wHat ##

    # Compute the integral along I_j of f(y_i | s; theta)
    Xb <- X%*%betaHat
    fval <- integration(fun = function(s, i,...) dfam(y[i], Xb[i] + gammaHat*s, phiHat, famly, size),
                        intmap, Kappas)
    if(!is.na(fval$warn) & PRINT!=0) cat(fval$warn, fill=T)

    # Self-consistent equations
    lwmat <- matrix(rep(Ipoint + (1-Ipoint)*(intmap[2,]-intmap[1,]), n), byrow = TRUE, nrow = n)
    dmat <- fval$results
    repeat{
      wOld <- wHat
      wmat <- matrix(rep(wHat, n), byrow = TRUE, nrow = n)
      numerator <- Kappas * dmat * wmat / lwmat
      denominator <- matrix(rep(rowSums(numerator), m), nrow=n)
      nuumat <- numerator/denominator
      colnames(nuumat) <- colnames(numerator)
      wHat <- colSums(nuumat, na.rm = TRUE)/n
      if (sum((wHat - wOld)^2) / sum(wOld^2) < tolpar)
        break
    }


    ## Update thetHat ##

    # Likelihood
    likfun <- function(x){
      lx <- length(x)
      if(PRINT!=0) cat("Likfun for x=",x, fill=T)

      if(x[lx] <= 0) x[lx] = 1e-18
      Xb <- X%*%x[1:lbet]
      eta <- function(i,s){Xb[i] + x[lx-1]*s}

      fval <- integration(fun = function(s, i,...) dfam(y[i], eta(i,s), x[lx], famly, size),
                          intmap, Kappas)
      if(!is.na(fval$warn) & PRINT!=0) cat(fval$warn, fill=T)

      o <- - sum(log(rowSums(  (Kappas * fval$results * wmat) / lwmat  ) +1e-18), na.rm=TRUE)
      if(PRINT!=0) cat("----> o=",o, fill=T)

      return(o)
    }

    likfun_s <- function(x){
      lx <- length(x)
      if(x[lx] <= 0) x[lx] = 1e-18
      Xb <- X%*%x[1:lbet]
      eta <- function(i,s){Xb[i] + x[lx-1]*s}

      fval <- integration(fun = function(s, i,...) dfam(y[i], eta(i,s), x[lx], famly, size),
                          intmap, Kappas, rel.tol = 1e-5)
      if(!is.na(fval$warn) & PRINT!=0) cat(fval$warn, fill=T)

      o <- - sum(log(rowSums(  (Kappas * fval$results * wmat) / lwmat  ) +1e-18), na.rm=TRUE)

      return(o)
    }

    grfun <- function(x){
      lx <- length(x)
      if(PRINT!=0) cat(" grfun for x=",x, fill=T)
      if(x[lx] <= 0) x[lx] = 1e-18

      Xb <- X%*%x[1:lbet]
      eta <- function(i,s){Xb[i] + x[lx-1]*s}
      mu <- function(i,s){gr_linkinv(eta(i,s))}
      mu.eta <- function(i,s){gr_mu.eta(eta(i,s))}


      fungrval1 <- function(i,s){
        ifelse(dfam(y[i], mu=mu(i,s), phi=x[lx], famly=famly, size=size)==0, 0,
               (((y[i]-mu(i,s)) * mu.eta(i,s)) / famly$variance(mu(i,s))) * 1 * dfam(y[i], mu=mu(i,s), phi=x[lx], famly=famly, size=size)
        )}
      fungrval2 <- function(i,s){
        ifelse(dfam(y[i], mu=mu(i,s), phi=x[lx], famly=famly, size=size)==0, 0,
               (((y[i]-mu(i,s)) * mu.eta(i,s)) / famly$variance(mu(i,s))) * s * dfam(y[i], mu=mu(i,s), phi=x[lx], famly=famly, size=size)
        )}


      fval <- integration(fun = function(s, i,...) dfam(y[i], eta(i,s), x[lx], famly, size),
                          intmap, Kappas, rel.tol = 1e-8, relTol = 1e-8)
      if(!is.na(fval$warn) & PRINT!=0) cat(fval$warn, fill=T)

      grval1 <- integration(fun = function(s, i,...) fungrval1(i,s),
                            intmap, Kappas, rel.tol = 1e-8, relTol = 1e-8)
      if(!is.na(grval1$warn) & PRINT!=0) cat(grval1$warn, fill=T)

      grval2 <- integration(fun = function(s, i,...) fungrval2(i,s),
                            intmap, Kappas, rel.tol = 1e-8, relTol = 1e-8)
      if(!is.na(grval2$warn) & PRINT!=0) cat(grval2$warn, fill=T)


      # numerator
      gr1_shared <- Kappas * grval1$results * wmat / lwmat
      gr1_terms <- sapply(1:ncol(X), function(j) {
        rowSums(gr1_shared * matrix(rep(X[,j],m), ncol=m))
      })
      gr2_term <- rowSums((Kappas * grval2$results * wmat) / lwmat)
      numerator <- cbind(gr1_terms, gr2_term)

      # denominator
      denominator <- rowSums(Kappas * fval$results * wmat / lwmat)
      denominator <- do.call(cbind, rep(list(denominator), ncol(numerator)))
      g <- -1/x[lx] * colSums(ifelse(denominator!=0,numerator/denominator,0), na.rm = TRUE)
      g <- c(g, nloptr::nl.grad(x[lx], function(phi) likfun_s(c(x[-lx], phi))))
      if(PRINT!=0) cat(" ----> g=",g, fill=T)
      if(PRINT==2){cat(" ----> numeric=",numDeriv::grad(likfun_s, x), fill=T)}

      return(g)
    }



    ### LBFGS given analytical gradient
    lowerin <- setNames(rep(-Inf, lbet+2), names(thetHat))
    lowerin["disp"] = 1e-18
    bbmle::parnames(likfun) <- names(thetHat)
    bbmle::parnames(grfun) <- names(thetHat)

    res <- bbmle::mle2(minuslogl = likfun,
                       gr = grfun,
                       start = thetHat,
                       lower = lowerin,
                       vecpar = TRUE,
                       method = "L-BFGS-B",
                       trace = (PRINT!=0))

    likv <- - res@details$value
    thetHat <- res@details$par
    betaHat <- thetHat[1:lbet]
    gammaHat <- thetHat[lbet+1]
    phiHat <- thetHat[lbet+2]

    # Relative difference with previous iteration
    reldif_par <- sum((wHat - wOld.out)^2) / sum(wOld.out^2) + sum((thetHat-thetOld.out)^2) / sum(thetOld.out^2)
    reldif_lik <- (likv - likvOld)^2/likvOld^2
    if(PRINT!=0) cat(paste0("Relative diffs: pars=", reldif_par,"  -  lik=", reldif_lik), fill=T)

    # save
    lreldif_par <- c(lreldif_par, reldif_par)
    lreldif_lik <- c(lreldif_lik, reldif_lik)
    lthetHat <- rbind(lthetHat, thetHat)
    llikv <- c(llikv, likv)
    cnv <- c(cnv, res@details$convergence)
    cnvmss <- c(cnvmss, res@details$message)

    # Stopping criteria
    if (reldif_par < tolpar) break
    if (reldif_lik < tolik){
      itolik = 1
      break
    }
    if (iter > maxiter){
      imaxi = 1
      break
    }
    if(difftime(Sys.time(), gotime, units = "sec") > maxsec){
      imaxt = 1
      break
    }
  }



  # Return
  cbi <- cbind(Tj = names(wHat), qj = intmap[1,], pj = intmap[2,], dens = wHat, distr=cumsum(wHat))
  rownames(cbi) <- NULL

  temp <- cbind(lthetHat,llikv)
  colnames(temp) <- c(colnames(lthetHat),"likv")

  tempr <- cbind(lreldif_par, lreldif_lik)
  colnames(tempr) <- c("par","likv")

  def_args <- formals(gelc.fit)
  call_args <- match.call(expand.dots=T)[-1]
  for (nm in names(call_args)){
    if(nm %in% c("y","zl","zr","X")) next
    if(nm == "famly"){
      def_args[[nm]] <- famly[1]
      next
    }
    def_args[[nm]] <- eval(call_args[[nm]], envir = parent.frame())
  }
  full_call <- as.call(c(as.name("gelc.fit"), def_args))

  return(list(thetaHat = thetHat, wHat = cbi, A = Kappas, likval = likv,
              convOPT = cnv, mssOPT = cnvmss,
              reachmaxi = imaxi, reachmaxt = imaxt, reachtolik = itolik,
              iter_vals = temp, iter_reldif = tempr,
              vcov = res@vcov,
              start_used = start_used, call = full_call))
}
