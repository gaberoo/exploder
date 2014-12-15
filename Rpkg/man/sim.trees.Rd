\name{sim.trees}
\alias{sim.trees}
\title{
Simulates sampled transmission trees.
}
\description{
Simulate an instance of a sampled transmission tree from epidemiological
parameters.
}
\usage{
sim.trees(N, beta, mu, psi, max.samples, min.outbreak = min(10,max.samples),max.time = -1.0) 
}
\arguments{
  \item{N}{
    Population size (carrying capacity). Set N=0 for N=Inf.
  }
  \item{beta}{
    Transmission rate.
  }
  \item{mu}{
    Death rate.
  }
  \item{psi}{
    Sampling rate.
  }
  \item{max.samples}{
    Simulate epidemic until this many individuals are sampled. Set to zero for
    species phylogenies (see Details).
  }
  \item{min.outbreak}{
    Minimum number of infecteds that is considered an outbreak (avoids
    stochastic extinctions).
  }
  \item{max.time}{
    Maximum time to run epidemic, even if max.samples isn't reached.
  }
}
\details{
  Fixed-time sampling.
  The default behavior is to simulate sampled phylogenies for sequentially
  sampled data. It is possible to simulate fixed time sampling by setting
  psi=0. In this case the simulation is run until the total number of
  infecteds is equal to max.samples, when the simulation is stopped and all
  lineages are sampled.
}
\value{
  Data frame with the event times and event types (0 = sampling, 1 =
  branching) in columns.
}
\references{
  Leventhal, Bonhoeffer & Stadler, 2012
}
\author{
  Gabriel E Leventhal
}
\seealso{
  \code{\link{expoTree}}
}
\examples{
  N <- 500
  beta <- 1
  mu <- 0.08
  psi <- 0.02
  nsamp <- 100

  tree <- sim.trees(N,beta,mu,psi,nsamp)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
