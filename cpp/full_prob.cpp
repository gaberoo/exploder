/*
 * Copyright (c) 2012-2014, Gabriel Leventhal, ETH Zurich
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

#include "expoTree.h"
#include "expoTreePSave.h"

#include "Forest.h"
#include "ExpoTreePars.h"

#include <gsl/gsl_statistics.h>

void printHelp();

// ===========================================================================

vector<double> K;        /* total population size */
vector<double> beta;     /* infection rate */
vector<double> mu;       /* recovery rate */
vector<double> psi;      /* sampling rate */
vector<double> rShift;   /* time of rate shift */
double rho = 0.0;        /* initial sampling probability */
int root = 0;            /* add a root to the tree */
int n = 0;               /* number of times */
int vflag = 0;           /* vebose flag */
int rescale = 1;
int est_norm = 1;
int extant = 0;
int maxExtant = 0;
int maxEvals = 1000;
double fixedRatio = -1.0;
char jointType = 'x';
char logFun = 'D';
int max_shifts;
int survival = 0;
int nroot = 0;
int model_type = 1;

// ===========================================================================

int main(int argc, char* argv[]) {
  string fn;

  // =========================================================================

  int c;
  opterr = 0;
  while ((c = getopt(argc,argv,"n:N:l:vb:u:s:o:hr:R:J:S:LE:t:g:")) != -1) {
    switch (c) {
      case 'n': nroot = atoi(optarg); break;
      case 'N': K.push_back(atof(optarg)); break;
      case 'v': ++vflag; break;
      case 'b': beta.push_back(atof(optarg)); break;
      case 'u': mu.push_back(atof(optarg)); break;
      case 's': psi.push_back(atof(optarg)); break;
      case 'o': rho = atof(optarg); break;
      case 'r': rescale = atoi(optarg); break;
      case 'R': fixedRatio = atof(optarg); break;
      // case 'h': printHelp(); return 0;
      case 'J': jointType = optarg[0]; break;
      case 'S': rShift.push_back(atof(optarg)); break;
      case 'L': logFun = 'L'; break;
      case 'E': est_norm = atoi(optarg); break;
      case 't': model_type = atoi(optarg); break;
      case 'g': survival = atoi(optarg); break;
    }
  }

  // =========================================================================

  if (vflag) fprintf(stderr,"Reading times...\n");

  if (optind == argc) {
    fprintf(stderr,"Please provide times file.\n");
    return 0;
  }

  Forest T(argc-optind,argv+optind);
  for (size_t i(0); i < rShift.size(); ++i) T.addRateShift(rShift[i]);

  if (vflag) {
    printf("Extant species at time zero = %d\n",T.at(0)->extant);
    printf("Max extant species = %d\n",T.maxExtant);
  }

  if (fixedRatio > 0.0) {
    mu.resize(psi.size());
    for (size_t i(0); i < psi.size(); ++i) {
      mu[i] = psi[i]*(1.0/fixedRatio-1.0);
    }
  }

  /*
  for (size_t i(0); i < K.size(); ++i) {
    printf("%6.2f ",K.at(i));
    printf("%6.2f ",beta.at(i));
    printf("%6.2f ",mu.at(i));
    printf("%6.2f | ",psi.at(i));
    if (i > 0) printf("%6.2f\n",rShift.at(i-1));
    else printf("%6s\n","0");
  }
  */

  vector<double> pars(5);
  pars[0] = beta[0];
  pars[1] = mu[0];
  pars[2] = psi[0];
  pars[3] = rho;
  pars[4] = K[0];

  double fx = 0.0;
  double surv = 0.0;
  double torig = 0.0;

  double* psave = NULL;

  extant = T.at(0)->extant;
  maxExtant = T.at(0)->maxExtant;

  int np = T[0]->times.size();
  int ldp = ceil(K[0]+1);
  psave = new double[np*ldp];

  fx = expoTreePSave(K.size(),K.data(),beta.data(),mu.data(),psi.data(),
                     rho,T[0]->times,T[0]->ttypes,&np,&ldp,psave,
                     extant,est_norm,vflag,rescale,nroot);

  cerr << "Log-likelihood = " << fx << endl;

  double ll = -INFINITY;
  double maxll;
  for (int j = 1; j < ldp; ++j) {
    ll = log(1.0);
    printf("%f %d %d %f\n",T[0]->times[0],0,j,ll);
  }
  for (int i = 1; i < np; ++i) {
    maxll = gsl_stats_max(psave+(i-1)*ldp+1,1,ldp-1);
    for (int j = 1; j < ldp; ++j) {
      // ll = log(psave[i*ldp+j])+psave[i*ldp];
      ll = log(psave[(i-1)*ldp+j])-log(maxll);
      printf("%f %d %d %f\n",T[0]->times[i],i,j,ll);
    }
    printf("\n");
    fflush(stdout);
  }

  delete[] psave;

  return 0;
}

// ===========================================================================

