runExpoTreeExt <- function(K,lambda,mu,psi,rho,times,ttypes,
                           survival=TRUE,shifts=NULL,vflag=0,
                           return.full=FALSE,rescale=TRUE,
                           root.lineages=0,estimate.norm=FALSE)
{
  p <- .Call("expoTreeEvalExt",
             K=K,lambda=lambda,mu=mu,psi=psi,rho=rho,
             times=as.numeric(times),
             ttypes=as.integer(ttypes),
             add_par=as.integer(c(0,vflag,rescale,
                                  root.lineages,
                                  estimate.norm)))

  lik <- NULL
  if (return.full) { lik <- p } else { lik <- p[root.lineages+2] }
  return(lik)
}


