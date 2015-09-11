#### functions for calculating local constant HDS(t) and its standard errors
#### functions for calculating beta(t), Lambda(t), and beta(t) SE in cox_lc.R
#### local constant HDS(t) SE basically just depends on beta(t) SE

#' Hazard discrimination summary estimate (local constant) at one time point
#'
#' \code{hdslc.fast} estimates HDS at a single time using the local-in-time
#' proportional hazards model. See Cai and Sun (2003, Scandinavian Journal of
#' Statistics) for details on the local-in-time PH model.
#'
#' The user typically will not interact with this function. Rather, \code{hdslc}
#' wraps this function and is what the user typically will use.
#'
#' @param S A vector of length \code{nrow(m)} (which is typically the number of
#' observations n), where each value is the subject-specific survival at time t
#' where t is implied by the choice of \code{betahat}.
#' @param betahat A p x 1 vector of coefficient estimates at time t of interest
#'   from the local-in-time Cox model. Vector length p is the number of
#'   covariates. Typically the output from \code{hdslc::finda} is passed here.
#' @param m A numeric n x p matrix of covariate values, with a column for each
#'   covariate and each observation is on a separate row.
#' @return The HDS estimate at times t, where t is implied by choice of \code{S}
#' and \code{betahat} passed to \code{hdslc.fast}.
hdslc.fast <- function(S, betahat, m){
  if(is.null(dim(m))) m <- matrix(m)
  else                m <- matrix(m, ncol=ncol(m))
  bm     <- m%*%betahat
  S3hat  <- sum(exp(2*bm)*S)
  S1hat  <- sum(S)
  S2hat  <- sum(exp(bm)*S)
  return(S3hat*S1hat/(S2hat^2))
}

#' Hazard discrimination summary (local constant) standard error estimate
#'
#' \code{hdslcse.fast} calculates an estimate of the variance for the
#' local constant hazard discrimination summary estimator at a time t. The time
#' t is implied by \code{S}, \code{betahat}, and \code{betahatse}
#'
#' The use will typically not interact with this function directly. Instead this
#' function is wrapped by \code{hdslc}.
#'
#' @param S A vector of length \code{nrow(m)} (which is typically the number of
#' observations n), where each value is the subject-specific survival at time t
#' where t is implied by the choice of \code{betahat}.
#' @param betahat A p x 1 vector of coefficient estimates at time t of interest
#'   from the local-in-time Cox model. Vector length p is the number of
#'   covariates. Typically the output from \code{hdslc::finda} is passed here.
#' @param m A numeric n x p matrix of covariate values, with a column for each
#'   covariate and each observation is on a separate row.
#' @param betahatse A p x p covariance matrix for betahat at time t
#' @return Variance estimate that has not been normalized. To get a usable
#'   standard error estimate, divide the output of this function by the
#'   bandwidth and sample size, and then take the square root.
hdslcse.fast <- function(S, betahat, m, betahatse){
  if(is.null(dim(m))) m <- matrix(m)
  else                m <- matrix(m, ncol=ncol(m))
  n    <- nrow(m)
  bm   <- m %*% betahat
  bmm  <- bm %*% rep(1, ncol(m))
  Smm  <- S %o% rep(1, ncol(m))

  S2hat <- colSums(2*m*exp(2*bmm)*Smm)/n
  S1hat <- colSums(m*exp(bmm)*Smm)/n
  #S0hat <- colSums(Smm)/n
  S0hat <- rep(0, ncol(m))

  Smat  <- rbind(S0hat, S1hat, S2hat)
  Sigma2 <- Smat %*% betahatse %*% t(Smat)

  Pf2  <- sum(exp(2*bm)*S)/n
  Pf0  <- sum(S)/n
  Pf1  <- sum(exp(bm)*S)/n

  A <- matrix(c(Pf2/(Pf1^2), -(2*Pf2*Pf0)/(Pf1^3), Pf0/(Pf1^2)), nrow=1)
  return(A %*% Sigma2 %*% t(A))
}

#' Hazard discrimination summary estimator
#'
#' \code{hdslc} returns local constant HDS estimates at all specified evaluation
#' times
#'
#' A local constant version of \code{hds}. While \code{hds} estimates HDS(t)
#' assuming the Cox proportional hazards model, \code{hdslc} estimates HDS(t)
#' using a relaxed, local-in-time Cox model. Specifically, the hazard
#' ratios are allowed to vary with time. See Cai and Sun (2003, Scandinavian
#' Journal of Statistics) and Tian Zucker Wei (2005, JASA) for details on the
#' local-in-time Cox model.
#'
#' Point estimates use \code{hdslc.fast} and standard errors use
#' \code{hdslcse.fast}. \code{hdslc.fast} requires an estimate of beta(t) (time-varying
#' hazard ratio), which is estimated using \code{finda()}; and subject specific
#' survival, which is estimated using sssf.fast(). \code{hdslcse.fast} the same
#' and in addition standard error estimates of beta(t), which are estimated
#' using \code{betahatse.fast()}.
#'
#' The covariate values \code{m} are centered for numerical stability. This is
#' particularly relevant for the standard error calculations.
#'
#' @param times A vector of observed follow up times.
#' @param status A vector of status indicators, usually 0=alive, 1=dead.
#' @param m A matrix or data frame of covariate values, with a column for each
#'   covariate and each observation is on a separate row. Non-numeric values
#'   are acceptable, as the values will be transformed into a numeric model
#'   matrix through \code{survival::coxph}.
#' @param evaltimes A vector of times at which to estimate HDS. Defaults to all
#'   the times specified by the \code{times} vector. If there are a lot of
#'   observations, then you may want to enter in a sparser vector of times for
#'   faster computation.
#' @param h A single numeric value representing the bandwdith to use, on the
#'   time scale. The default bandwidth is a very ad hoc estimate using
#'   Silverman's rule of thumb
#' @param se TRUE or FALSE. TRUE: calculate and return standard error estimates.
#'   FALSE: do not calculate standard errors estimates and return NAs. Defaults
#'   to TRUE. May want to set to FALSE to save computation time if using this
#'   function to compute bootstrap standard errors.
#' @examples
#' hdslc(times = survival::pbc[1:312, 2],
#'       status = (survival::pbc[1:312, 3]==2)*1,
#'       m = survival::pbc[1:312, 5])
#' @return A data frame with three columns: 1) the evaluation times, 2) the HDS
#'   estimates at each evaluation time, and 3) the standard error estimates at
#'   each evaluation time
#' @export
#' @importFrom survival Surv coxph
hdslc <- function(times,
                  status,
                  m,
                  evaltimes=times[order(times)],
                  h=1.06*sd(times)*(length(times)^(-0.2)),
                  se=TRUE){
  m   <- as.matrix(m)
  fit <- coxph(Surv(times, status)~m, x = TRUE)
  m   <- fit$x
  m   <- scale(m, scale = FALSE)     # center marker values
  m   <- as.matrix(m)

  timesort   <- order(times)     # this and next three lines sort data by time
  m_s        <- m[timesort, , drop=FALSE]    # sorted marker values
  status_s   <- status[timesort] # sorted status values
  times_s    <- times[timesort]  # sorted times

  evalt <- findInterval(evaltimes, times_s) # find indices of evaluation times
  evaln <- length(evalt)                    # how many evaluation times

  betaD <- apply(matrix(times_s[evalt]), 1, function(x)
    finda(x, times=times_s, status=status_s, covars=m_s, h=h))
  betaD <- t(betaD)
  if(nrow(betaD)==1) betaD <- t(betaD)
  es <- sssf.fast(betat = betaD, status = status_s, m = m_s, evalt)
  hdslcres <- apply(matrix(1:evaln), 1, function(x)
    hdslc.fast(es[x, ], betaD[x, ], m_s))
  hdslcseres <- rep(NA, evaln)
  if(se){
    # calculate standard errors if se=TRUE
    betase <- betahatse.fast(betahat = betaD, times = times_s, status = status_s,
                             m = m_s, h = h, evalt)
    hdslcseres <- apply(matrix(1:evaln), 1, function(x){
      hdslcse.fast(es[x,], betaD[x, ], m_s, betase[x,,])})
    hdslcseres <- sqrt(hdslcseres/(h*length(times_s)))
  }
  return(data.frame(times=evaltimes, hdslchat=hdslcres, se=hdslcseres))
}

