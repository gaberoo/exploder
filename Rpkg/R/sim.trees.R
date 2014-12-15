sim.trees <- function(N,beta,mu,psi,max.samples,
                      min.outbreak=min(10,max.samples),max.time=-1.0) 
{
  len <- 2*max.samples
  times <- ttypes <- vector(length=len)
  out <- .C("RSimTrees",
            N=as.integer(N),
            beta=as.numeric(beta),
            mu=as.numeric(mu),
            psi=as.numeric(psi),
            max_samples=as.integer(max.samples),
            min_outbreak=as.integer(min.outbreak),
            max_time=as.numeric(max.time),
            maxlen=as.integer(len),
            times=as.numeric(times),
            ttypes=as.integer(ttypes))
  out <- data.frame(times=out$times,ttypes=out$ttypes)
  out <- out[order(out$times),]
  out$lineages <- cumsum(1-2*out$ttypes)
  out
}

