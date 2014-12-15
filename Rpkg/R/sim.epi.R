sim.epi <- function(N,beta,mu,psi,max.samples,
                    min.outbreak=min(10,max.samples),
                    model.type=1,
                    max.time=-1.0,
                    old.state=list()) 
{
  pars <- c(N,beta,mu,psi)
  out <- .Call("RSimEpi",parameters=pars,
               max_samples=max.samples,
               min_outbreak=min.outbreak,
               model_type=model.type,
               max_time=max.time,state=old.state)
  o <- order(out$times)
  out$times <- out$times[o]
  out$ttypes <- out$ttypes[o]
  out
}


