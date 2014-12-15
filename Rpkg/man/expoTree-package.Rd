\name{expoTree-package}
\alias{expoTree-package}
\alias{expoTree}
\docType{package}
\title{
  Calculate density dependent likelihood of a phylogenetic tree
}
\description{
  Calculates the density dependent likelihood of a phylogenetic tree. 
  It takes branching and sampling times as an argument and integrates the
  likelihoood function over the whole tree.
}
\details{
  \tabular{ll}{
    Package: \tab expoTree\cr
    Type: \tab Package\cr
    Version: \tab 1.0.1\cr
    Date: \tab 2013-09-03\cr
    License: \tab BSD 3 clause\cr
  }
}
\author{
  Gabriel E Leventhal, partly adapted from MATLAB code by Awad H. Al-Mohy and
  using the routines DLNAC1 and DLARPC by Sheung Hun Cheng, and DLAPST from
  ScaLAPACK.
}
\references{
  Leventhal, Guenthard, Bonhoeffer & Stadler, (2013) "Using an epidemiological
  model for phylogenetic inference reveals density-dependence in HIV
  transmission".

  Al-Mohy & Higham (2011) "Computing the Action of the Matrix Exponential, with
  an Application to Exponential Integrators". SIAM Journal on Scientific
  Computing, 33.

  Cheng & Highami (2001) "Implementation for LAPACK of a Block Algorithm for
  Matrix 1-Norm Estimation". Numerical Analysis Report No. 393, Manchester
  Centre for Computational Mathematics, Manchester, England, August 2001, and
  LAPACK Working Note 152.

  Blackford et al. (1997) "ScaLAPACK Users' Guide". Society for Industrial and
  Applied Mathematics, Philadelphia, PA.
}

\keyword{phylodynamics}
\examples{
  # simulate trees
  N <- 15
  beta <- 1
  mu <- 0.1
  psi <- 0.1
  rho <- 0
  nsamp <- 20
  forest <- lapply(1:2,function(x) sim.trees(N,beta,mu,psi,nsamp))

  # plot lineages-through-time
  plotLTT(forest)

  extant <- sapply(forest,function(t) sum(2*t[,2]-1))
  lineages <- lapply(forest,function(t) sum(2*t[,2]-1)+cumsum(1-2*t[,2]))
  max.lineages <- sapply(lineages,max)

  # calculate likelihood for the forest
  lik <- sapply(forest,function(tree) {
    runExpoTree(matrix(c(N,beta,mu,psi,rho),nrow=1),tree[,1],tree[,2])
  })
  cat("Likelihood = ",sum(lik),"\n")

  if (! is.nan(lik)) {
    # Find optimal parameters
    opt <- expoTree.optim(forest,
                   lo=c(max(max.lineages),1e-10,1e-10),hi=c(20,3,1),
                   fix=c(FALSE,FALSE,TRUE,FALSE,TRUE),
                   fix.val=c(0,0,-psi/(psi+mu),0,0),
                   control=list(trace=1,REPORT=1,maxit=10))
  }
}
