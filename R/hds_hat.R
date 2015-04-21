#' Hazard discrimination summary estimate at one time point
#'
#' \code{hds_t} estimates HDS at time t under the PH assumption
#'
#' The user typically will not interact with this function. Rather, \code{hds}
#' is a wrapper for this function and is what the user typically will use.
#'
#' @param t The time at which to calculate HDS
#' @param L0hat A data frame with variable names of hazard and time. Typically
#'   the object returned by \code{basehaz}.
#' @param betahat A vector of coefficient estimates from the Cox model.
#'   Typically the \code{coefficients} value from the \code{coxph.object} object
#'   returned by \code{coxph}.
#' @param m A numeric matrix of covariate values, with a column for each
#'   covariate and each observation is on a separate row.
#' @importFrom survival Surv coxph basehaz
hds_t <- function(t, L0hat, betahat, m){
  if(is.null(dim(m))) m <- matrix(m)
  else                m <- matrix(m, ncol=ncol(m))
  L      <- L0hat$hazard[findInterval(t, L0hat$time)]
  bm     <- m%*%betahat
  S3hat  <- sum(exp(2*bm - exp(bm)*L))
  S1hat  <- sum(exp(-exp(bm)*L))
  S2hat  <- sum(exp(bm - exp(bm)*L))
  return(S3hat*S1hat/(S2hat^2))
}


#' Hazard discrimination summary estimator
#'
#' \code{hds} returns HDS estimates at all specified evaluation times
#'
#' A wrapper for \code{hds_t}. Since \code{hds_t} only estimates HDS at one time
#' point, this function calls \code{hds_t} multiple times to estimate the entire
#' HDS curve. It is the main function the user will interact with in this
#' package.
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
#' @examples
#' hds(times = survival::pbc[1:312, 2], status = (survival::pbc[1:312, 3]==2)*1,
#'     m = survival::pbc[1:312, 5])
#' @return A data frame with three columns: 1) the evaluation times, 2) the HDS
#'   estimates at each evaluation time, and 3) the standard error estimates at
#'   each evaluation time
#' @export
#' @importFrom survival Surv coxph basehaz
hds <- function(times, status, m, evaltimes=times[order(times)]){
#  for(i in 1:ncol(m)){
    # center the covariates
#    m[,i] <- m[,i] - mean(m[,i])
#  }
  fit <- coxph(Surv(times, status)~m, x = TRUE)
  L0  <- basehaz(fit)
  m <- fit$x
  m <- scale(m, scale = FALSE)     # center marker values
  timesort   <- order(times)     # this and next three lines sort data by time
  m_s        <- m[timesort, ]    # sorted marker values
  status_s   <- status[timesort] # sorted status values
  times_s    <- times[timesort]  # sorted times
#  fit        <- coxph(Surv(times_s, status_s)~m_s) # fit PH model
#  L0         <- basehaz(fit, centered=TRUE)  # calculate baseline cumulative hazard
  pt_est <- apply(matrix(evaltimes), 1,
                 function(x) hds_t(t=x, L0hat=L0, betahat=fit$coef, m=m_s)) # calculate HDS(t) at all times
  pt_se  <- hds_se(time=times_s, status=status_s, m=m_s, evaltimes=evaltimes)
  return(data.frame(times=evaltimes, hdshat=pt_est, se=pt_se))
}


