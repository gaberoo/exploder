/*
 * Copyright (c) 2014, Gabriel Leventhal, ETH Zurich
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of the ETH Zurich nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
using namespace std;

#include "array.h"
#include "expoTreeExt.h"
#include "Forest.h"
#include "ExpoTreePars.h"

// ===========================================================================

int root = 0;            /* add a root to the tree */
int n = 0;               /* number of times */
int vflag = 2;           /* vebose flag */
int rescale = 1;
int est_norm = 1;
int extant = 0;
int maxExtant = 0;
int maxEvals = 1000;
double fixedRatio = -1.0;
char jointType = 'g';
char logFun = 'D';
int max_shifts;
int survival = 0;
int nroot = 0;

int main(int argc, char* argv[]) {
  if (vflag) fprintf(stderr,"Reading times...\n");

  if (optind == argc) {
    fprintf(stderr,"Please provide times file.\n");
    return 0;
  }

  Forest* T = new Forest(argc-optind,argv+optind);

  if (vflag) {
    printf("Extant species at time zero = %d\n",T->at(0)->extant);
    printf("Max extant species = %d\n",T->maxExtant);
  }

  int N = 100;
  vector<double> K(1,1.0*N);
  Array2D<double> lambda(1,N+1);
  Array2D<double> mu(1,N+1);
  Array2D<double> psi(1,N+1);

  for (int i = 0; i <= N; ++i) {
    lambda[i] = 1.0;
    mu[i] = 0.1*i/N;
    psi[i] = 0.1;
  }

  double rho = 0.0;

  double fx = expoTreeEvalExt(K,lambda,mu,psi,rho, 
                              (*T)[0]->times,(*T)[0]->ttypes,
                              (*T)[0]->extant+nroot,
                              est_norm,vflag,rescale,
                              nroot);

  cout << fx << endl;

  delete T;
  return 0;
}
