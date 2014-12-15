lik.bdsky <- function(pars,times,ttypes,shifts,survival=TRUE,nroot=0,vflag=0,
                      correct.ordering=FALSE) {
  .Call("bdSkyEval",parameters=pars,times=times,ttypes=ttypes,shifts=shifts,
        add_par=as.integer(c(survival,vflag,nroot)))
}

lik.bdsky.2 <- function(pars,times,ttypes,shifts,survival=TRUE,nroot=0,vflag=0,
                        correct.ordering=FALSE) {
# sanity checks
  m <- length(shifts)+1
  if (m < nrow(pars)) { return(-Inf) }
# make sure times are sorted
  ot <- order(times)
  times <- times[ot]
  ttypes <- ttypes[ot]
  shifts <- sort(shifts)
# get parameters ready
  x <- cbind(pars,c(shifts,max(times)))
# calculate lineages through time at time shifts
  lineages <- .get.lineages(cbind(times,ttypes),nroot=nroot)
  shift.lin <- approx(times,lineages,x[,5],"constant",rule=1,f=1)$y
# get intervals of all times
  ival <- findInterval(times,shifts)+1
# recursively calculate Bi and pi(t[i-1])
  ABp <- .bdsky.ABp(x,m,vflag)
  Ai <- ABp$Ai
  Bi <- ABp$Bi
  pi.ti <- ABp$pi.ti
# calculate likelihood for the branching and sampling events
  lik <- vector("numeric",length(times))
  ln.Qi <- vector("numeric",length(times))
  for (i in 1:length(times)) {
    k <- ival[i]
    ti <- ifelse(k>1,-x[k-1,5],0.0)
    ln.Qi[i] <- .bdsky.lnq(t0=-times[i],ti=ti,Ai=Ai[k],Bi=Bi[k])
    if (vflag > 0) { 
      message(sprintf("ti(%d) = % f, lnq = % f",k,ti,ln.Qi[i])) 
    }
    if (ttypes[i] == 1) {
      lik[i] <- log(2) + log(x[k,1]) + ln.Qi[i]  # infections
    } else if (ttypes[i] == 0) {
      lik[i] <- log(x[k,3]) - ln.Qi[i]  # samplings
    }
  }
  log.lik <- sum(lik)
  if (vflag > 0) { message(sprintf("loglik = %f\n",log.lik)) }
# sampling + shifts
  ti <- 0
  ln.X <- vector("numeric",m)
  for (k in 1:(m-1)) {  ########## <- why don't we include the last branch ?
  #for (k in 1:m) {
    t0 <- -x[k,5]
    lnq <- .bdsky.lnq(t0,ti,Ai[k],Bi[k])
    ln.X[k] <- shift.lin[k]*lnq
    if (vflag > 0) {
      message(sprintf("lnq(%.1f|%d) = % f, lins = %d",t0,k,lnq,shift.lin[k]))
    }
    ti <- t0
  }
  if (vflag) {
    message(sprintf("log-lik[trans] = %f",sum(lik[ttypes==1])))
    message(sprintf("log-lik[sampl] = %f",sum(lik[ttypes==0])))
    message(sprintf("log-lik[shift] = %f",sum(ln.X)))
  }
  log.lik <- log.lik + sum(ln.X)
  log.lik <- log.lik - log(2*x[m,1]) ########## <- I don't understand this correction
  if (vflag > 0) { message(sprintf("loglik = %f\n",log.lik)) }
  if (survival) {
    t0 <- -x[m,5]
    ti <- ifelse(m>1,-x[m-1,5],0.0)
    ll.surv <- log(1-.bdsky.p(x[m,],t0,ti,Ai[m],Bi[m]))
    if (vflag) {
      message(sprintf("log-lik[surv ] = %f",ll.surv))
    }
    log.lik <- log.lik - ll.surv
  }
  if (correct.ordering) {
    log.lik <- log.lik + log(2)*(sum(ttypes==1)-1)
  }
  unname(log.lik)
}

.bdsky.ABp <- function(x,m,vflag=0) {
  Ai <- apply(x,1,.bdsky.A)
  Bi <- vector("numeric",m)
  pi.ti <- vector("numeric",m+1)
  ti <- 0         # ti = tm
  pi.ti[1] <- 1   # p[m+1](tm)
  for (i in 1:m) {
    t0 <- -x[i,5]
    Bi[i] <- .bdsky.B(x[i,],Ai[i],pi.ti[i])
    pi.ti[i+1] <- .bdsky.p(x[i,],t0,ti,Ai[i],Bi[i])
    if (vflag > 0) {
      message(sprintf("A(%d) = % f, B(%d) = % f, p(%d) = % f",
                      i,Ai[i],i,Bi[i],i+1,pi.ti[i+1]))
    }
    ti <- t0
  }
  list(Ai=Ai,Bi=Bi,pi.ti=pi.ti)
}

.bdsky.A <- function(x) {
  sqrt((x[1]-x[2]-x[3])**2 + 4*x[1]*x[3])
}

.bdsky.B <- function(x,Ai,pi.neg1) {
  ((1-2*(1-x[4])*pi.neg1)*x[1]+x[2]+x[3])/Ai
}

.bdsky.p <- function(x,t0,ti,Ai,Bi) {
  a <- x[1]+x[2]+x[3]
  eA <- exp(Ai*(ti-t0))
  num <- eA*(1+Bi)-(1-Bi)
  den <- eA*(1+Bi)+(1-Bi)
  (a-Ai*num/den)/(2*x[1])
}

.bdsky.q <- function(t0,ti,Ai,Bi) {
  eA <- exp(-Ai*(t0-ti))
  4*eA/((eA*(1+Bi)+(1-Bi))**2)
}

# ti is the upper bound
.bdsky.lnq <- function(t0,ti,Ai,Bi) {
  log.eA <- Ai*(ti-t0)
  log(4.0) + log.eA - 2*(log.eA + log(1+Bi) + log( 1+(1-Bi)/(exp(log.eA)*(1+Bi)) ))
}

