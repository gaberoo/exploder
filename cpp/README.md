# expoTree: calculate the likelihood of a phylogenetic tree with density-dependent transmission

This is the C++ implementation of the code used to calculate the likelihood of
a phylogenetic tree given an epidemiological susceptible-infected (SI) transmission model.

Cite this as Leventhal, Bonhoeffer, Günthardt and Stadler, 2013 (doi:XXXXXXX)

This code is distributed according to the BSD Licence.

## Dependencies

A C/C++ and FORTRAN compiler are required to compile the code.

### FORTRAN compiler

On Mac OS X, you can get one from one of these sources:

* [HPC Mac OS X](http://hpc.sourceforge.net/)
* [R for Mac OS X](http://cran.r-project.org/bin/macosx/tools/)

### Library dependencies:

* BLAS (distributed with Mac OS X via the Accelerate Framework)
* GNU Scientific Library (GSL). You can get it here: http://www.gnu.org/software/gsl/
* MPI for compiling parallel versions of the code

## Compiling instructions

First, on the command line, clone the repository.

    :::bash
    git clone https://bitbucket.org/gaberoo/sitree.git

Edit the file "Make.inc" such that the environment variables point to the
C/C++ and FORTRAN compilers on your system. An example is provided in the file
"Make.inc.example". This should work for Mac OS X systems, provided a FORTRAN
compiler is installed.

Then compile the source.

    :::bash
    make

## Running instructions

A phylogenetic tree is represented by it's branching and sampling times.
These are provided to the programs in text files (one per tree), with the
event time in the first column and the event type (0 = sampling, 1 =
branching) in the second column. The event times do not need to be in order,
as there's sorted when read.

To calculate the likelihood of a set of trees (a forest) using the
density-dependent SI model with parameters "N", "beta", "mu" and "psi", run,

    :::bash
    ./calc_likelihood -N <populationSize> -b <beta> -u <mu> -s <psi> tree1.txt [tree2.txt] [tree3.txt] ...

To run an MCMC integration using the DREAM algorithm, run,

    :::bash
    ./dream -N <maxN> -b <betaMin> -B <betaMax> -u <muMin> -U <muMax> \
       -s <psiMin> -S <psiMax> -o 0 -e <maxEvals> -t <lockCode> \
       -G <gelmanEvals> -X <burnIn> -J j -g <gelmanCrit> \
       -C <outputFilePrefix> tree1.txt [tree2.txt] [tree3.txt] ...

* `xxxMin` and `xxxMax` are the lower and upper limits of the prior 
  distribution of parameters. The lower limit of `N` is determined from the 
  maximum number of extant lineages at any timepoint in the tree.
* `maxEvals` are the maximum number of likelihood evaluations to perform. 
* The `lockCode` is a binary code that locks certain parameters to a value:
    * N = 16
    * beta = 8
    * mu = 4
    * psi = 2
    * rho = 1
  For example, if you want to lock beta and psi: lockCode = 8 + 2 = 10.  
  **For epidemiological inference, "rho" should always be locked to rho = 0, 
  using `-o 0 -t 1` or another lock code!**
* `gelmanEvals` sets the number of function evaluations between which
  convergence is evaluated using the Gelman-Rubin criterion.
* `burnIn` are the number of function evaluations during burn-in. Note that in
  the DREAM algorithm calibration is performed during burn-in, such that
  detailed balance is not guaranteed.
* `gelmanCrit` is the critical value for the potential scale reduction factor.
  Default is 1.05.
* `outputFilePrefix` is the prefix for the output files. Each of the 5 chains
  will be stored in a file outputFilePrefix-X.txt

## Source files

* *Makefile, Make.inc*:
  Makefiles to compile code. Edit *Make.inc* to fit your system. The
  provided file *should* work on Mac OS X, provided gcc, g++ and gfortran are
  available.

* *shared*:
  This directory contains a C implementation of the Al-Mohy & Higham MATLAB
  code from 2010. Uses the DLACN1 algorithm written in FORTRAN 66.
  **Reference:** A. H. Al-Mohy and N. J. Higham, Computing the action of
  the matrix exponential, with an application to exponential
  integrators. MIMS EPrint 2010.30, The University of Manchester, 2010.

* *expomv.c*:
  Performs the likelihood calculation described in Leventhal, Bonhoeffer,
  Günthardt and Stadler, 2013 (in prep) for an array of branching/sampling
  times.

* *expoTree.h, expoTree.cpp*:
  C++ helper functions to calculate the likelihood of a phylogenetic tree given an SI
  model.
  
* *Tree.h, Forest.h*:
  C++ classes for the storage of phylogenetic trees represented by their
  branching and sampling times.

* *calc_likelihood.cpp*:
  Calculate the likelihood of a tree given epidemiological parameters and
  branching/sampling times.

* *dream.cpp*:
  Runs a Bayesian MCMC simulation using the DREAM algorithm
  **Reference**: Vrugt, J. A., ter Braak, C. J. F., Diks, C. G. H., Robinson, B. A., Hyman, 
  J. M., Higdon, D., 2009. Accelerating Markov chain Monte Carlo simulation by differential 
  evolution with self-adaptive randomized subspace sampling. International Journal of Nonlinear 
  Sciences and Numerical Simulation 10 (3), 273-290. doi:10.1515/IJNSNS.2009.10.3.273

* *sim_trees.cpp*:
  Simulate an SI epidemic for given parameters *beta, mu, psi* and *N*. Can
  output the phylogenetic tree or the epidemic incidence.

## Event types

0. sampling-removal event -- apply function `expo.fs`
1. branching event -- apply function `expo.ft`
2. sampling-only event -- apply function `expo.fs2`
4. discovery event -- increases lineage count, but does not apply function
20. rate shift event
21. rate shift event
22. rate shift event

## Licence

Copyright (c) 2011-2013 Institute of Integrative Biology, ETH Zurich. 
                        All rights reserved.

Copyright (c) 2010      Nick Higham and Awad Al-Mohy

Copyright (c) 1992-2011 The University of Tennessee and The University
                        of Tennessee Research Foundation.  All rights
                        reserved.

Copyright (c) 2000-2011 The University of California Berkeley. All
                        rights reserved.

Copyright (c) 2006-2011 The University of Colorado Denver.  All rights
                        reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

The copyright holders provide no reassurances that the source code
provided does not infringe any patent, copyright, or any other
intellectual property rights of third parties.  The copyright holders
disclaim any liability to any recipient for claims brought against
recipient by any third party for infringement of that parties
intellectual property rights.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

