## examples for local-in-time Cox model

#' @importFrom stats runif
rt <- function(x1,x2){
  # generates random survival times using Cai and Sun (2003) simulation example
  (1/x1)*log(-(x1*log(runif(length(x1)))/(exp(0.5*x2)))+1)
}

# true HDS under Cai and Sun (2003) data generating example where
# lambda(t|x1,x2) = exp(t*x1 + 0.5*x2)
lmbd <- function(t, x1, x2){
  # conditional hazard function
  exp(t*x1 + 0.5*x2)
}

S <- function(t, x1, x2){
  # conditional survival function
  exp(-exp(0.5*x2) * ((exp(t*x1)-1)/x1))
}

# hdstruecs <- rep(NA, length(seq(0,2,0.01)))
# i <- 0
# for(tt in seq(0, 2, 0.01)){
#   H2 <- 0
#   H1 <- 0
#   H0 <- 0
#   i <- i+1
#   for(x1 in seq(-1,1,length.out=500)){
#     for(x2 in seq(-4,4,length.out=500)){
#       l <- lmbd(tt, x1, x2)
#       s <- S(tt, x1, x2)
#       fx2 <- dnorm(x2)
#       H2 <- H2 + (l*l*s*fx2)
#       H1 <- H1 + (l*s*fx2)
#       H0 <- H0 + (s*fx2)
#     }
#   }
#   hdstruecs[i] <- H2 * H0 / H1^2
#   cat(tt, " ", hdstruecs[i], '\n')
# }

#' @importFrom stats runif rnorm
do.one.cs.fast <- function(n=500, h=0.5, evaltimes="fast"){
  # evaltimes="fast" defaults to calculating at only 100 evenly spaced
  # times between 0.02 and 2
  # otherwise put in your own vector of evaluation times (e.g. all observed times)
  x1 <- runif(n, 1, 3)
  x2 <- rnorm(n)
  D <- cbind(rt(x1,x2), rep(1, length(x1)), x1, x2)
  D <- D[which(!is.na(D[,1])), ]
  D <- D[order(D[,1]),]
  if(evaltimes == "fast"){
    evalt <- findInterval(seq(0.02, 2, 0.02), D[,1])
  } else{
    evalt <- findInterval(evaltimes, D[,1])
  }
  evaln <- length(evalt)
  betaD <- apply(matrix(D[evalt,1]), 1, function(x) finda(x, times=D[,1], status=D[,2], covars=D[,3:4], start=c(0,0), h=h))
  betaD <- t(betaD)
  cat("beta done; ")
  es <- sssf.fast(betat = betaD, status = D[,2], m = D[,3:4], evalt)
  cat("S done; ")
  hdslcres <- apply(matrix(1:evaln), 1, function(x) hdslc.fast(es[x, ], betaD[x, ], D[,3:4]))
  cat("hds done; ")
  betase <- betahatse.fast(betahat = betaD, times = D[,1], status = D[,2], m = D[,3:4], h = h, evalt)
  cat("beta SE done; ")
  hdslcseres <- apply(matrix(1:evaln), 1, function(x){
    hdslcse.fast(es[x,], betaD[x, ], D[,3:4], betase[x,,])})
  cat("hds SE done", "\n")
  #  par(mfrow=c(3,1))
  #  plot(seq(0,2,0.01), hdstruecs, ylim=c(0,2.5), xlim=c(0,2), type="n")
  #  lines(seq(0,2,0.01), hdstruecs, ylim=c(0,2.5), xlim=c(0,2), lwd=2)
  #  lines(D[,1], hdslcres, col='red')
  #  lines(D[,1], hdslcres+1.96*sqrt(hdslcseres/(h*nrow(D))), lty=2, col='red')
  #  lines(D[,1], hdslcres-1.96*sqrt(hdslcseres/(h*nrow(D))), lty=2, col='red')
  #  plot(D[,1], betaD[,1], ylim=c(-0.5,1.5), xlim=c(0,2), type='n')
  #  lines(D[,1], betaD[,1], col='red')
  #  abline(0,1, lwd=2)
  #  lines(D[,1], betaD[,1]+1.96*sqrt(betase[,1,1]/(h*nrow(D))), col='red', lty=2)
  #  lines(D[,1], betaD[,1]-1.96*sqrt(betase[,1,1]/(h*nrow(D))), col='red', lty=2)
  #  plot(D[,1], betaD[,2], ylim=c(0,1), xlim=c(0,2), type='n')
  #  lines(D[,1], betaD[,2], col='red')
  #  abline(h=0.5, lwd=2)
  #  lines(D[,1], betaD[,2]+1.96*sqrt(betase[,2,2]/(h*nrow(D))), col='red', lty=2)
  #  lines(D[,1], betaD[,2]-1.96*sqrt(betase[,2,2]/(h*nrow(D))), col='red', lty=2)
  #  dev.off()
  return(rbind(D[evalt,1], hdslcres, hdslcseres, t(betaD), betase[,1,1], betase[,2,2]))
}

L0.fast <- function(betat, status, m, evalt){
  # calculates just the cumulative baseline hazard
  # this is just used for model checking
  n <- length(status)
  evaln <- length(evalt)
  #evalt <- round(seq(1, n, length.out=100))
  etbm <- exp(tcrossprod(betat, m))
  etbmn <- etbm
  doit <- rep(NA, evaln)
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
  cumsum(doit)
}
