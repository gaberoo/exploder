\name{plotLTT}
\alias{plotLTT}
\title{
Plot lineages through time.
}
\description{
Creates a (or adds to an existing) plot of the lineages-through-time of a
list of phylogenetic trees.
}
\usage{
plotLTT(trees, col=rgb(.5,.5,.5,.5), xlab="Time", ylab="Lineages", log="y", add=FALSE, 
        plot = TRUE, plot.tree = NA, ...) 
}
\arguments{
  \item{trees}{
    A list of trees in two-column format. Column 1 is the branching/sampling
    time and column two is the event type (1 = branching, 0 = sampling) 
  }
  \item{col}{
    Color of the lineages.
  }
  \item{xlab}{
    x-axis label.
  }
  \item{ylab}{
    y-axis label.
  }
  \item{log}{
    Use logarithmic axis (defaults to logarithmic y-axis).
  }
  \item{add}{
    Add the lineages to the current plotting device.
  }
  \item{plot}{
    If set to FALSE, only the number of lineages at each time point are
    returned.
  }
  \item{plot.tree}{
    Plot the phylogenetic tree in the background. Only works for one tree.
  }
  \item{...}{
    Other parameters to pass to 'plot()'.
  }
}
\references{
  Leventhal, Bonhoeffer, Guenthard & Stadler, 2012
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

  trees <- lapply(1:10,function(i) sim.trees(N,beta,mu,psi,nsamp))
  plotLTT(trees)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
