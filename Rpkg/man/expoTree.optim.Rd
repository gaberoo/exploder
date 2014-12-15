\name{expoTree.optim}
\alias{expoTree.optim}
\title{
Perform optimization to find a maximum likelihood estimate.
}
\description{
Perform optimization to find a maximum likelihood estimate.
}
\usage{
expoTree.optim(forest, lo = rep(5,0), hi = rep(5,0),
      fix = rep(5,F), fix.val = rep(5,0), pars = vector(length = sum(!fix))*1,
      survival = TRUE, method = "pso", control = list(trace = 0))
}
\arguments{
  \item{forest}{
    List of trees in two-column format. Column 1 are the event times
    (branching or sampling) and column 2 is the event type (0 = sampling, 1 =
    branching).
  }
  \item{lo}{
    Lower bound for parameters.
  }
  \item{hi}{
    Upper bound for parameters.
  }
  \item{fix}{
    Logical vector specifying which parameters to keep constant.
  }
  \item{fix.val}{
    Values for fixed parameters. Also specify dummy values for variable
    parameters. These are ignored.
  }
  \item{pars}{
    Starting values for the parameters (only required for Nelder-Mead
    optimization).
  }
  \item{survival}{
    Condition on the likelihood of observing the tree.
  }
  \item{method}{
    Choose optimization method:
    pso         = use Particle Swarm optimization (R 'pso' package);
    nelder-mead = use the Nelder-Mead algorithm from the 'optim' method. 
  }
  \item{control}{
    Control parameters for the optimization method. See 'psoptim' and 'optim' for details.
  }
}
\details{
  Special case:
  It is possible to fix the ratio between mu and psi. If mu is fixed at a
  negtive value, it is seen as a fixed ratio.
  Example: 
  r <- 0.5
  fix <- c(F,F,T,F,T)
  fix.val <- c(0,0,-r,0,0)
  At every evaluation, mu is calculated accordingly:
  psi/(psi+mu) = r  ==> mu = psi*(1/r-1)
}
\value{
  Output from the optimization. List with the following entries:
    par        : parameter estimates
    value      : negative log-likelihood
  For the other returned values, see help(psoptim) or help(optim) respectively
}
\references{
  Leventhal, Bonhoeffer, Guenthard & Stadler, 2012
}
\author{
  Gabriel E. Leventhal
}
\seealso{
  \code{\link{expoTree}}
}
\examples{
  # simulate trees
  N <- 15
  beta <- 1
  mu <- 0.1
  psi <- 0.1
  rho <- 0
  nsamp <- 20
  tree <- sim.trees(N,beta,mu,psi,nsamp)

  extant <- sum(2*tree[,2]-1)
  lineages <- sum(2*tree[,2]-1)+cumsum(1-2*tree[,2])
  max.lineages <- max(lineages)

  # calculate likelihood for the forest
  lik <- runExpoTree(matrix(c(N,beta,mu,psi,rho),nrow=1),tree[,1],tree[,2])
  cat("Likelihood = ",sum(lik),"\n")

  if (! is.nan(lik)) {
    expoTree.optim(list(tree),
                   lo=c(max(max.lineages),1e-10,1e-10),hi=c(50,3,1),
                   fix=c(FALSE,FALSE,TRUE,FALSE,TRUE),
                   fix.val=c(0,0,-psi/(psi+mu),0,0),
                   control=list(trace=1,REPORT=1,maxit=10))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{likelihood}
\keyword{phylogenetics}
