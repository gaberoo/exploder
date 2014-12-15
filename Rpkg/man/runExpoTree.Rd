\name{runExpoTree}
\alias{runExpoTree}
\title{
Density-dependent likelihood caluculation
}
\description{
Calculates the density dependent likelihood of a phylogenetic tree. 
It takes branching and sampling times as an argument and integrates the
likelihoood function over the whole tree.
}
\usage{
runExpoTree(pars, times, ttypes, survival = TRUE, shifts = NULL, vflag = 0)
}
\arguments{
  \item{pars}{
    Matrix of parameters. If the number of columns is larger than 5, only the
    first 5 columns are used. Each row is of the form c(N,beta,mu,psi,rho).
    The rows correspond to rate shifts.
  }
  \item{times}{
    Vector of event times.
  }
  \item{ttypes}{
    Vector of event types.
  }
  \item{survival}{
    Condition on the likelihood of observing the tree.
  }
  \item{shifts}{
    Times of rate shifts.
  }
  \item{vflag}{
    Verbosity level.
  }
}
\value{
  Log-Likelihood
}
\details{
  Parameter matrix:

  N = Total population size. Must be at least the maximal number of extant lineages
  at any given time.

  beta = Branching rate. For epidemic trees, this is the infection rate. For
  species trees, this is the speciation rate.

  mu = Extinction rate. For epidemic trees, this is the death/recovery rate.
  For species trees, this is the extinction rate.

  psi = Sampling rate.

  rho = Initial sampling rate.
}
\references{
  Leventhal, Guenthard, Bonhoeffer & Stadler, 2013
}
\author{
  Gabriel E Leventhal
}
\seealso{
  \code{\link{expoTree}}
}
\examples{
times <- cumsum(runif(10))
ttypes <- rep(1,10)

N <- 15
beta <- 1
mu <- 0.1
psi <- 0
rho <- 1
pars <- matrix(c(N,beta,mu,psi,rho),nrow=1)
lik <- runExpoTree(pars,times,ttypes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
