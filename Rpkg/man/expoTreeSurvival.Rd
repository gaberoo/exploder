\name{expoTreeSurvival}
\alias{expoTreeSurvival}
\title{
  Calculate the probability that at least one lineage is sampled.
}
\description{
  Calculate the probability that at least one lineage is sampled.
}
\usage{
expoTreeSurvival(pars, times, ttypes, shifts = NULL, vflag = 0) 
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
  \item{shifts}{
    Times of rate shifts.
  }
  \item{vflag}{
    Verbosity level.
  }
}
\value{
  Likelihood
}
\references{
  Tanja Stadler, Journal of Theoretical Biology, 2010
}
\author{
  Gabriel E Leventhal
}

\seealso{
\code{\link{expoTree}}
}
\examples{
  N <- 15
  beta <- 1
  mu <- 0.1
  psi <- 0
  rho <- 1
  torig <- 10.0

  pars <- matrix(c(N,beta,mu,psi,rho),nrow=1)

  lik <- expoTreeSurvival(pars,times=torig,ttypes=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
