runExpoTree <- function(pars,times,ttypes,survival=TRUE,shifts=NULL,vflag=0,
                        return.full=FALSE,rescale=TRUE,root.lineages=0,
                        root.correction=FALSE,estimate.norm=FALSE)
{
  if (! is.matrix(pars)) pars <- matrix(as.numeric(pars),nrow=1)
  matDim <- dim(pars)
  nshifts <- length(shifts) + sum(ttypes %in% 20:22)

  if (matDim[2] < 5) {
    cat("Minimum 5 columns required in pars.\n")
    return(-Inf)
  }

  if (matDim[1] < nshifts+1) {
    cat("Not enough parameters for shifts.\n")
    return(-Inf)
  }

  if (matDim[1] > nshifts+1) pars <- pars[1:(nshifts+1),]

  if (! is.null(shifts)) {
    times <- c(times,shifts)
    ttypes <- c(ttypes,rep(20,length(shifts)))
    o <- order(times)
    times <- times[o]
    ttypes <- ttypes[o]
  }

  p <- .Call("expoTreeEval",
             parameters=pars,
             times=as.numeric(times),
             ttypes=as.integer(ttypes),
             add_par=as.integer(c(0,vflag,rescale,
                                  root.lineages,
                                  root.correction,
                                  estimate.norm)))
  # p = (p0,p1,p2,...,pN)
  if (return.full) {
    lik <- p
  } else {
    lik <- p[root.lineages+2]
  }

  if (survival) {
    surv <- expoTreeSurvival(pars,times,ttypes,vflag=vflag,
                             return.full=return.full,rescale=rescale,
                             estimate.norm=TRUE)
    if (return.full) surv <- surv[2]
    if (surv > -Inf) {
      lik <- lik - surv
    } else {
      lik <- -Inf
    }
  }

  return(lik)
}


