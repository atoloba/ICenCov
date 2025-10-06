#' Construct augmented Turnbull intervals and incidence matrix 
#'
#' @param L,R Numeric vectors of left/right endpoints of the observed intervals.
#' @param Lin,Rin Logical scalar or length-\eqn{n} vectors indicating whether
#'   the left/right endpoints are included. 
#' @param xmin,xmax Numeric bounds of the support of \eqn{Z}. 
#'   (defaults: \code{xmin = 0}, \code{xmax = max(R)}).
#' @param keepall Logical; if \code{FALSE} (default), remove partition intervals
#'   with zero contribution (columns with all zeros). 
#'
#' @details
#' This function adapts function Aintmap from the \CRANpkg{interval} package 
#' to handle augmented Turnbull intervals. 
#'
#' @seealso \code{\link[interval]{Aintmap}}
#' @export
Aintmap2 <- function(L,R, Lin=TRUE, Rin=TRUE, xmin=0, xmax=max(R), keepall=FALSE){
  n <- length(L)
  if(length(Lin)==1) Lin <- rep(Lin,n)
  if(length(Rin)==1) Rin <- rep(Rin,n)
  
  # Treat [L,R) as (L-eps,R), and (L,R] as (L,R+eps)
  LRvalues <- sort(unique(c(xmin, L, R, xmax)))
  eps <- min(diff(LRvalues)[diff(LRvalues) > .Machine$double.eps^0.5]) / 3
  Le <- L
  Re <- R
  Le[Lin] <- L[Lin] - eps
  Re[Rin] <- R[Rin] + eps
  oLR <- order(c(Le, Re+eps/2))
  cLR <- c(L,R)[oLR]
  cLRin <- c(Lin,Rin)[oLR]
  
  # label L=1 and R=0
  Leq1.Req0 <- c(rep(1,n), rep(0,n))
  cLeq1.Req0 <- Leq1.Req0[oLR]
  
  # Delete observations that are entering or leaving the interval at the exact same time.
  rmove <- duplicated(cbind(cLR, cLeq1.Req0, cLRin))
  cLR <- cLR[!rmove]
  cLeq1.Req0 <- cLeq1.Req0[!rmove]
  cLRin <- cLRin[!rmove]
  
  # Delete those overlapping L-open with R-closed, and L-closed with R-open
  rmove <- ifelse(cLR[-length(cLR)] == cLR[-1] & 
                    cLeq1.Req0[-length(cLR)] == !cLeq1.Req0[-1] & 
                    cLRin[-length(cLR)] == !cLRin[-1], T, F)
  cLR <- cLR[!rmove]
  cLeq1.Req0 <- cLeq1.Req0[!rmove]
  cLRin <- cLRin[!rmove] 
  
  # Construction of intervals: e.g. case of Lin=T, Rin=T 
  # intmapL is closed when an observation is entering the interval and open if it's leaving,
  # intmapR is the other way around.
  # The idea is to keep the observation inside intmapL, intmapR
  intmapL <- cLR[-length(cLR)]
  intmapR <- cLR[-1]
  intmapLin <- vector(length = length(cLR)-1)
  intmapRin <- vector(length = length(cLR)-1)
  intmapLin[as.logical(cLeq1.Req0[-length(cLR)])] <- cLRin[-length(cLR)][as.logical(cLeq1.Req0[-length(cLR)])]
  intmapLin[!as.logical(cLeq1.Req0[-length(cLR)])] <- !cLRin[-length(cLR)][!as.logical(cLeq1.Req0[-length(cLR)])]
  intmapRin[as.logical(cLeq1.Req0[-1])] <- !cLRin[-1][as.logical(cLeq1.Req0[-1])]
  intmapRin[!as.logical(cLeq1.Req0[-1])] <- cLRin[-1][!as.logical(cLeq1.Req0[-1])]
  intmapview <- cbind(intmapL, intmapR, intmapLin, intmapRin)
  # View(intmapview)
  
  # verification
  check <- intmapview[intmapL==intmapR & (intmapLin!=1 | intmapRin!=1)]
  if(length(check)>0) stop("intmap nondefined")
  
  intmap <- matrix(c(intmapL, intmapR), byrow=TRUE, nrow=2)
  attr(intmap,"LRin") <- matrix(c(intmapLin, intmapRin), byrow=TRUE, nrow=2)
  m <- dim(intmap)[[2]]
  Lbracket <- rep("(", m)
  Lbracket[intmapLin] <- "["
  Rbracket <- rep(")", m)
  Rbracket[intmapRin] <- "]"
  intname <- paste(Lbracket,intmapL,",",intmapR,Rbracket,sep="")
  
  # Incidence matrix
  Le[!Lin] <- L[!Lin] - eps
  Re[!Rin] <- R[!Rin] + eps
  A <- matrix(0,n,m,dimnames=list(1:n,intname))
  for (i in 1:n){
    intmapLe <- intmapL
    intmapRe <- intmapR
    if(Lin[i]==T) intmapLe[intmapLin] <- intmapLe[intmapLin]+eps
    if(Rin[i]==T) intmapRe[intmapRin] <- intmapRe[intmapRin]-eps
    tempint <- Le[i]<intmapLe & intmapRe<Re[i]
    A[i,tempint]<-1
  }
  
  if(!keepall){
  # Reduction of Kappa dimension: some intervals have no contribution from observations, so w = 0
  idel <- which(colSums(A) == 0)
  if(length(idel)>0){
    A <- A[,-idel]
    intmap <- intmap[,-idel]
    attr(intmap,"LRin") <- attr(intmap,"LRin")[,-idel]
  } 
  }
  
  if (any(rowSums(A)==0)) stop("Some individuals are not contributing to any Turnbull interval")
  return(mget(c("A", "intmap")))
}





#' Turnbull's estimator along the augmented Turnbull intervals 
#'
#' @details
#' This function is adapted from the \CRANpkg{interval} package to operate on
#' augmented Turnbull partitions. 
#'
#' @seealso \code{\link[interval]{icfit}}
#' @export
icfit2 <- function(L, R, Lin=T, Rin=T, keepall=FALSE, accu=1e-10){
  Aintmap <- Aintmap2(L,R,Lin,Rin,keepall=keepall)
  A <- Aintmap$A
  n <- length(L)
  m <- ncol(A)
  wHat <- rep(1/m, m)
  iter = 0
  repeat{
    wOld <- wHat
    wmat <- matrix(rep(wHat, n), byrow = TRUE, nrow = n)
    numerator <- A * wmat
    denominator <- matrix(rep(rowSums(numerator), m), nrow=n)
    nuumat <- ifelse(denominator==0, 0, numerator/denominator)
    colnames(nuumat) <- colnames(numerator)
    wHat <- colSums(nuumat)/n
    iter = iter+1
    if (sum((wHat - wOld)^2) / sum(wOld^2) < accu)
      break
  }
  intmap <- Aintmap$intmap
  pf <- wHat
  return(mget(c("intmap","pf","A")))
}