\name{infExpoTree}
\alias{infExpoTree}
\title{
Likelihood calculation without density dependence.
}
\description{
Calculates the likelihood of a phylogenetic tree without density dependence. 
It takes branching and sampling times as an argument and integrates the
likelihoood function over the whole tree.
}
\usage{
infExpoTree(beta, mu, psi, rho, times, ttypes, survival = TRUE)
}
\arguments{
  \item{beta}{
    Branching rate. For epidemic trees, this is the infection rate. For species
    trees, this is the speciation rate.
  }
  \item{mu}{
    Extinction rate. For epidemic trees, this is the death/recovery rate. For species
    trees, this is the extinction rate.
  }
  \item{psi}{
    Sampling rate.
  }
  \item{rho}{
    Initial sampling rate.
  }
  \item{times}{
    Vector of event times.
  }
  \item{ttypes}{
    Vector of event types.
  }
  \item{survival}{
    Condition on survival of the tree.
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
  times <- cumsum(runif(10))
  ttypes <- rep(1,10)

  N <- 15
  beta <- 1
  mu <- 0.1
  psi <- 0
  rho <- 1

  lik <- infExpoTree(beta,mu,psi,rho,times,ttypes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
