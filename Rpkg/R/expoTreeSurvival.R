expoTreeSurvival <- function(pars,times,ttypes,shifts=NULL,vflag=0,
                             return.full=FALSE,rescale=TRUE,
                             root.lineages=0,estimate.norm=FALSE) 
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

  ttypes[! ttypes %in% 20:22] <- 99

  if (! is.null(shifts)) {
    times <- c(times,shifts)
    ttypes <- c(ttypes,rep(20,length(shifts)))
    o <- order(times)
    times <- times[o]
    ttypes <- ttypes[o]
  }

  p <- .Call("expoTreeEval",parameters=pars,
             times=as.numeric(times),
             ttypes=as.integer(ttypes),
             survival=as.integer(c(1,vflag,rescale,root.lineages,0,
                                   estimate.norm)))

  # p = (p0,p1,p2,...,pN)
  if (return.full) {
    return(p)
  } else {
    return(p[2])
  }
}

