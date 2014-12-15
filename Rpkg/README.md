# expoTree

expoTree is an R packge to calculate the likelihood of the a phylogenetic
tree, under the assumption of a stochastic susceptible-infected model for the
underlying population dynamics.

expoTree is not yet available on CRAN, so you have to download and
install it "manually". Make sure that you have a C/C++ and a FORTRAN compiler
installed.

## Installation

### CRAN installation

**expoTree** is now available on CRAN! Enter

    :::R
    install.packages("expoTree")

in an R console or choose from the list of available packages.

### Source installation

Download the latest release from the "Downloads" section or clone the Git 
repository to get the latest version.

    :::bash
    git clone https://bitbucket.org/gaberoo/expotree.git

If you downloaded a package, type the following on the command line

    :::bash
    R CMD INSTALL expoTree_{currentVersion}.tar.gz

Replace {currentVersion} with whatever version you installed. If you clone the Git repo,

    :::bash
    R CMD INSTALL expotree

You can now use the library in R.

## Usage

Check the R help files for a detailed description. Here's a quick-and-dirty
example:

    :::R
    library(expoTree)

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

## MCMC and expoTree

You can of course use expoTree together with MCMC packages in R. However, I
strongly suggest that you use the [pure C version](https://bitbucket.org/gaberoo/expotree) 
if you plan on performing Bayesian inference. It conatins a C implementation
of the DREAM algorithm 
doi:[10.1515/IJNSNS.2009.10.3.273](http://dx.doi.org/10.1515/IJNSNS.2009.10.3.273]).
The DREAM algorithm speeds of convergence of the Metropolis-Hastings
algorithm, which is great because the calculation of the likelihood is really
expensive.

