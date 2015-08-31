cVV  <- function(time, m, status, betahat, betavar, n){
  # used by the hds_se() function to calculate the covariance between beta and Lambda_0
  # returns an n x (p+1) matrix, where p is number of covariates
  if(is.null(dim(m))) m <- matrix(m, ncol=1)
  z0e <- exp(m%*%betahat) %*% rep(1, ncol(m))
  z1e <- m*z0e
  z2e <- m*z1e
  time <- matrix(time, ncol=1)
  status <- matrix(status*1, ncol=1)

  num <- matrix(NA, nrow=length(time), ncol=ncol(m))
  den <- num
  V1  <- matrix(NA, nrow=length(time), ncol=1)
  count <- 0
  for(i in time){
    count <- count+1
    timei <- (time>=i) %*% rep(1, ncol(m))
    num[count, ] <- n*colSums(z1e*timei)/colSums(z0e*timei)^2
    den[count, ] <- (colSums(z2e*timei)/colSums(z0e*timei)) - (colSums(z1e*timei)/colSums(z0e*timei))^2
    V1[count, 1]  <- n*n/sum(z0e[,1]*timei[,1])^2
  }
  num <- num[order(time),]
  num <- as.matrix(num)
  num[status[order(time)]==0, ] <- 0
  cnum <- apply(num, 2, cumsum)
  return(cbind((1/n)*cumsum(V1[order(time)]*status[order(time)]) + (1/n)*diag(cnum%*%betavar%*%t(cnum)),
               cnum%*%betavar))
}

hds_se <- function(time, status, m, evaltimes){
  fit     <- coxph(Surv(time, status)~m)
  betahat <- matrix(fit$coef, ncol=1)
  betavar <- fit$var
  L0hat   <- basehaz(fit, centered=TRUE)
  cVVout  <- cVV(time, m, status, betahat, betavar, fit$n) ### come back to this
  se      <- rep(NA, length(evaltimes))

  if(is.null(dim(m))) m <- matrix(m, ncol=1)
  m    <- matrix(m, ncol=ncol(m))
  bm   <- m%*%betahat
  bmm  <- bm %*% rep(1, ncol(m))
  ebmm <- exp(bmm)
  for(i in 1:length(evaltimes)){
    if(evaltimes[i]<min(L0hat$time) | evaltimes[i]>max(L0hat$time)){
      # if evaluation time is outside of observed times, set to NA
      se[i]    <- NA
    } else{
      Lh <- L0hat$hazard[findInterval(evaltimes[i], L0hat$time)]
      C  <- exp(2*bmm - ebmm*Lh)
      B  <- exp(bmm - ebmm*Lh)
      A  <- exp(-ebmm*Lh)
      Cb <- C*(2*m-(m*ebmm)*Lh) # n x p matrix, where p is number of covariates
      Bb <- B*(m-(m*ebmm)*Lh)   # n x p matrix, where p is number of covariates
      Ab <- A*((-m*ebmm)*Lh)    # n x p matrix, where p is number of covariates
      CL <- C*(-ebmm)
      BL <- B*(-ebmm)
      AL <- A*(-ebmm)

      dbeta   <- (colSums(Cb)*colSums(A)/colSums(B)^2) + (colSums(C)*colSums(Ab)/colSums(B)^2) -
        (2*colSums(C)*colSums(A)*colSums(Bb)/colSums(B)^3)
      # length p vector
      dlambda <- (colSums(CL)*colSums(A)/colSums(B)^2) + (colSums(C)*colSums(AL)/colSums(B)^2) -
        (2*colSums(C)*colSums(A)*colSums(BL)/colSums(B)^3)

      #    var1    <- dbeta%*%(fit$var*fit$n)%*%dbeta + dlambda[1]^2*(cVVout[i,1]) + sum(2*dbeta*dlambda*(cVVout[i,2:ncol(cVVout)]))
      cVVind  <- findInterval(evaltimes[i], time)
      var1    <- dbeta%*%(fit$var*fit$n)%*%dbeta + dlambda[1]^2*(cVVout[cVVind,1]) +
        sum(2*dbeta*dlambda*(cVVout[cVVind,2:ncol(cVVout)]))
      # we could do this with one big matrix multiplication as well
      #  var1    <- dbeta^2*(fit$var*fit$n) + dlambda^2*(cVVout[1,]) + 2*dbeta*dlambda*(cVVout[2,])

      m11 <- mean(A[,1]^2)-mean(A[,1])^2
      m12 <- mean(A[,1]*B[,1])-mean(A[,1])*mean(B[,1])
      m13 <- mean(A[,1]*C[,1])-mean(A[,1])*mean(C[,1])
      m22 <- mean(B[,1]^2)-mean(B[,1])^2
      m23 <- mean(B[,1]*C[,1])-mean(B[,1])*mean(C[,1])
      m33 <- mean(C[,1]^2)-mean(C[,1])^2

      sigma    <- matrix(c(m11, m12, m13, m12, m22, m23, m13, m23, m33), nrow=3)
      gradient <- c(mean(C[,1])/(mean(B[,1])^2), -2*mean(A[,1])*mean(C[,1])/(mean(B[,1])^3), mean(A[,1])/(mean(B[,1])^2))
      var2     <- gradient%*%sigma%*%gradient

      se[i]    <- sqrt((var1+var2)/fit$n)
      # we can do this b/c cov should be zero (asymptotically)
      #    cat(var1, var2, se[i], "\n")
    }
  }
  return(se)
}
