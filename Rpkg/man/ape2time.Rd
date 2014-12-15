\name{ape2time}
\alias{ape2time}
\title{
Convert tree from 'phylo' format (ape) to event time format.
}
\description{
Extract branching and sampling times from a tree in the native ape format.
}
\usage{
ape2time(tree, eps = 0.0)
}
\arguments{
  \item{tree}{
    Tree in 'phylo' format.
  }
  \item{eps}{
    Error tolerance for the present time (see Details).
  }
}
\details{
  In some cases R doesn't handle mathematical operations quite well and it can
  happen that tip ages are calculated as 10^{-14} or similar instead of 0.0.
  'eps' defines the tolerance for what is considered extant at the present
  (tip age < eps).
}
\author{
  Gabriel E Leventhal
}
\seealso{
  \code{\link{expoTree}}, \code{\link{ape}}
}
\examples{
  data(bird.families)
  ape2time(bird.families)
}
\keyword{phylogeny}
\keyword{phylo}
\keyword{ape}
