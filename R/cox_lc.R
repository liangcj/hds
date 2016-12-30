## collection of functions that implement estimation of beta(t) and subject
## specific survival (which requires estimation of cumulative baseline hazard)
## of local constant Cox model as specified by

## Cai and Sun (2003, Scan J of Stat) and Tian Zucker Wei (2005, JASA)
## finda(): local constant estimate of beta(t) using Cai and Sun (2003)
## betahatse.fast(): estimate SE of beta(t) using Cai and Sun (2003)
## sssf.fast(): subject specific survival at all m and t values
##              using Tian Zucker Wei (2005)

#' Estimate the time-varying coefficients from a local-in-time Cox model
#'
#' \code{finda} estimates the time-varying coefficients beta(t) at a single time
#' from a local-in-time Cox model. Think of it as a Cox model where the
#' the coefficients are allowed to vary with time. Further details can be found
#' in Cai and Sun (2003) and Tian et al. (2005).
#'
#' The naming of the function \code{finda} stands for "find a(t)", where "a(t)"
#' is the notation used in Cai and Sun (2003) to represent the time-varying
#' Cox model coefficients. We also refer to "a(t)" as "beta(t)" through the documentation.
#'
#' The user typically will not interact with this function, as \code{finda} is
#' wrapped by \code{hdslc}.
#'
#' @param tt Time to estimate beta(t) at
#' @param times A vector of observed follow up times.
#' @param status A vector of status indicators, usually 0=alive, 1=dead.
#' @param covars A matrix or data frame of numeric  covariate values, with a
#'   column for each covariate and each observation is on a separate row.
#' @param start A vector of length p of starting values to be passed to
#'   \code{optim} for the numerical optimization procedure. p is the number of
#'   covariates. Defaults to all zeroes.
#' @param h A single value on the time scale representing the bandwidth to use.
#' @param ... Additional parameters to pass to \code{optim}.
#' @importFrom stats optim
#' @return A vector of length p, where p is the number of covariates. The vector
#' is the estimated beta(t) from the local-in-time Cox model at time \code{tt}.
#' @references Cai Z and Sun Y (2003). Local linear estimation for time-dependent coefficients
#'   in Cox's regression models. \emph{Scandinavian Journal of Statistics}, 30: 93-111.
#'   \href{https://doi.org/10.1111/1467-9469.00320}{doi:10.1111/1467-9469.00320}
#' @references Tian L, Zucker D, and Wei LJ (2005). On the Cox model with time-varying
#'   regression coefficients. \emph{Journal of the American Statistical Association},
#'   100(469):172-83. \href{https://doi.org/10.1198/016214504000000845}{doi:10.1198/016214504000000845}
finda <- function(tt, times, status, covars, start=rep(0, ncol(covars)), h=400,...){
  # MAKE SURE DATA IS SORTED BY TIME
  # find a(t) for local constant Cox regression as per Cai and Sun (2003)
  # returns a(t) at t=tt
  # uses Epanechnikov kernel
  # added gradient function for roughly 10% speedup in numerical fitting
  nn   <- 1-((times-tt)/h)^2 >= 0 # negative values from epanechnikov
  K    <- (1-((times-tt)/h)^2)*0.75/h # epanechnikov kernel
  #  nn   <- (times-tt)>=(-h) & (times-tt)<=h # negative values from rectangular kernel
  #  K    <- 1/(2*h)                          # rectangular kernel
  N    <- length(times)
  llpl <- function(a){
    if(is.null(ncol(covars))) lp <- covars*a
    else           lp <- covars %*% a
    -sum((K*status*(lp - log(cumsum(exp(lp)[N:1])[N:1])))[nn])
  }
  llplg <- function(a){
    # gradient function
    if(is.null(ncol(covars))) lp <- covars*a
    else           lp <- covars %*% a
    elp <- exp(lp)
    num <- apply(sweep(covars, 1, elp, "*"), 2, function(x) rev(cumsum(rev(x))))
    den <- cumsum(elp[N:1])[N:1] %o% rep(1, ncol(covars))
    #-colSums((K*status*(covars - num/den))[nn,])
    #-colSums(sweep((covars - num/den), 1, K*status, "*")[nn,])
    -colSums(((covars - num/den) * (K*status))[nn, , drop=FALSE])
  }
  #  fit <- optim(par=start, fn=llpl, control=list(maxit=2000), ...)
  fit <- optim(par=start, fn=llpl, gr = llplg, method="BFGS", ...)
  #  cat(fit$convergence, sum(nn), "\n")
  return(fit$par)
}

#' Calculate estimates of the covariance matrix for beta(t).
#'
#' To help with overall computational efficiency of \code{hdslc}, this function
#' will return multiple covariance matrices - one covariance matrix for each
#' requested evaluation time.
#'
#' See Cai and Sun (2003, Scandinavian Journal of Statistics) for details on the
#' local constant estimator for beta(t). The user will not typically interact
#' with this function. The function is wrapped by \code{hdslc}, and the output
#' is used by \code{hdslcse.fast}.
#'
#' @return An N x p x p matrix, where N is the number of evaluation times and
#'   p is the number of covariates.
#' @keywords internal
#' @importFrom tensor tensor
betahatse.fast <- function(betahat, times, status, m, h, evalt){
  n  <- length(status)
  evaln <- length(evalt)
  #evalt <- round(seq(1, n, length.out=100))
  #betahat <- betahat[rep(1:100, c(diff(evalt), 1)), ]
  betahat <- betahat[rep(1:evaln, diff(c(0,evalt[1:(evaln-1)],n))), ]
  ebmm <- exp(tcrossprod(betahat, m))
  S0 <- apply(matrix(1:n), 1, function(x) sum(ebmm[x,x:n]))/n
  S1 <- matrix(NA, nrow=n, ncol=ncol(m))
  for(i in 1:n){S1[i,] <- (ebmm[i,i:n, drop=FALSE] %*% m[i:n, , drop=FALSE])/n}
  # can speed up by writing C code that does this
  mouter <- array(NA, dim=c(dim(m), ncol(m)))
  for(i in 1:n){mouter[i, , ] <- m[i, ] %o% m[i, ]}
  S2 <- array(NA, dim=c(n, ncol(m), ncol(m)))
  for(i in 1:n){S2[i, , ] <- tensor::tensor(ebmm[i,i:n, drop=FALSE], mouter[i:n, , , drop=FALSE], 2, 1)/n}
  A <- sweep(S2, 1, S0, "/")
  Btmp <- sweep(S1, 1, S0, "/")
  B <- array(NA, dim=c(dim(m), ncol(m)))
  for(i in 1:n){B[i, , ] <- outer(Btmp[i, ], Btmp[i, ])}
  V <- A-B
  tmat <- apply(matrix(times[evalt]), 1, function(x) times-x)
  # each i-th column is the n times minus the i-th time
  K <- ((1-((tmat)/h)^2)*0.75/h)*(1-((tmat)/h)^2 >= 0) # epanechnikov kernel
  integ <- tensor::tensor(V, K, 1, 1)/n
  I <- array(NA, dim=c(evaln, ncol(m), ncol(m)))
  for(i in 1:evaln){
    try(I[i, , ] <- solve(integ[ , , i]))
  }
  return(I * 0.6)
}

sssf.fast <- function(betat, status, m, evalt){
  # calculates subject specific survival for all t and all m values
  # faster, vectorized version
  n <- length(status)
  evaln <- length(evalt)
  #evalt <- round(seq(1, n, length.out=100))
  etbm <- exp(tcrossprod(betat, m))
  exp(-apply(matrix(1:n), 1, function(x){
    etbmn <- etbm / etbm[, x]
    #cumsum(status * (1/apply(matrix(1:n), 1, function(y) sum(etbmn[y, y:n]))))
    doit <- rep(NA, evaln)
    #for(i in 1:99){
    #  tmp <- sum(etbmn[i, (evalt[i]):n]) - c(0, cumsum(etbmn[i, (evalt[i]+1):(evalt[i+1]-1)]))
    #  doit[i] <- (1/tmp) %*% status[evalt[i]:(evalt[i+1]-1)]
    #}
    #doit[100] <- status[n] * (1/etbmn[i, n])
    for(i in 2:evaln){
      if(evalt[i]==0) doit[i] <- NA
      else{
        if(evalt[i]==1) doit[i] <- status[1] * (1/sum(etbmn[1, ]))
        else{
          if((evalt[i] == evalt[i-1]) | (evalt[i] == evalt[i-1]+1)) doit[i] <- (1/sum(etbmn[i, (evalt[i]):n])) * status[evalt[i]]
          else{
            tmp <- sum(etbmn[i, (evalt[i]):n]) + c(rev(cumsum(etbmn[i, (evalt[i]-1):(evalt[i-1]+1)])), 0)
            doit[i] <- (1/tmp) %*% status[(evalt[i-1]+1):(evalt[i])]
          }
        }
      }
    }
    if(evalt[1]==0) doit[1] <- NA
    else{
      if(evalt[1]==1) doit[1] <- (1/sum(etbmn[1, ])) * status[1]
      else{
        tmp <- sum(etbmn[1, (evalt[1]):n]) + c(rev(cumsum(etbmn[1, (evalt[1]-1):1])), 0)
        doit[1] <- (1/tmp) %*% status[1:(evalt[1])]
      }
    }
    #doit[1] <- (1/sum(etbmn[1, ])) * status[1]
    cumsum(doit)
  }))
}

