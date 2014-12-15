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

#include "expoTreeExt.h"
#include "Forest.h"
#include "ExpoTreePars.h"
#include "array.h"

#include "rapidjson/document.h"

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
char jointType = 'g';
char logFun = 'D';
int max_shifts;
int survival = 0;
int nroot = 0;

// ===========================================================================

int main(int argc, char* argv[]) {
  string user_fun;

  // =========================================================================

  int c;
  opterr = 0;
  while ((c = getopt(argc,argv,"n:N:l:vf:o:hr:R:S:LE:")) != -1) {
    switch (c) {
      case 'n': nroot = atoi(optarg); break;
      case 'N': K.push_back(atof(optarg)); break;
      case 'v': ++vflag; break;
      case 'f': user_fun = optarg; break;
      case 'o': rho = atof(optarg); break;
      case 'r': rescale = atoi(optarg); break;
      case 'h': printHelp(); return 0;
      case 'S': rShift.push_back(atof(optarg)); break;
      case 'L': logFun = 'L'; break;
      case 'E': est_norm = atoi(optarg); break;
    }
  }

  // =========================================================================

  if (vflag) fprintf(stderr,"Reading times...\n");

  if (optind == argc) {
    fprintf(stderr,"Please provide times file.\n");
    return 0;
  }

  Forest* T = new Forest(argc-optind,argv+optind);
  for (size_t i(0); i < rShift.size(); ++i) T->addRateShift(rShift[i]);

  if (vflag) {
    printf("Extant species at time zero = %d\n",T->at(0)->extant);
    printf("Max extant species = %d\n",T->maxExtant);
  }

  if (fixedRatio > 0.0) {
    mu.resize(psi.size());
    for (size_t i(0); i < psi.size(); ++i) {
      mu[i] = psi[i]*(1.0/fixedRatio-1.0);
    }
  }

  for (size_t i(0); i < K.size(); ++i) {
    printf("%6.2f ",K.at(i));
    if (i > 0) printf("%6.2f\n",rShift.at(i-1));
    else printf("%6s\n","0");
  }

  int N = 100;
  vector<double> K(1,1.0*N);
  Array2D<double> lambda(1,N+1);
  Array2D<double> mu(1,N+1);
  Array2D<double> psi(1,N+1);

  /* Reading JSON */
  rapidjson::Document d;
  d.Parse<0>(user_fun.c_str());

  rapidjson::SizeType i;
  rapidjson::Value& a = d["lambda"];
  assert(a.IsArray());
  for (i = 0; i < a.Size(); ++i) lambda[i] = a[i].GetDouble();

  a = d["mu"];
  assert(a.IsArray());
  for (i = 0; i < a.Size(); ++i) mu[i] = a[i].GetDouble();

  a = d["psi"];
  assert(a.IsArray());
  for (i = 0; i < a.Size(); ++i) psi[i] = a[i].GetDouble();

  double fx = expoTreeEvalExt(K,lambda,mu,psi,rho, 
                              (*T)[0]->times,(*T)[0]->ttypes,
                              (*T)[0]->extant+nroot,
                              est_norm,vflag,rescale,
                              nroot);

  cout << fx << endl;
  delete T;
  return 0;
}

// ===========================================================================

void printHelp() {
  printf("calc_likelihood: Calculate likelihood for a phylogenetic tree from an epidemic\n\n");
  printf("usage: calc_likelihood -N <population_size> -b <beta> -u <mu> -s <psi>\n");
  printf("                       [-rvh] <times_files> \n\n");
  printf("  N : population size\n");
  printf("  b : infection rate 'beta'\n");
  printf("  u : recovery rate 'mu'\n");
  printf("  s : sampling rate 'psi'\n");
  printf("  o : present sampling rate 'rho'\n");
  printf("  r : rescale staring vector after each iteration\n");
  printf("  l : model to use (0 = density independent; 1 = density dependent)\n");
  printf("  v : verbose (can be used multiple times)\n");
  printf("  h : print this help message\n");
  printf("  J : join type.\n");
  printf("        g (default) : joint likelihood\n");
  printf("        j           : joint likelihood + survival\n");
  printf("        m           : mean likelihood\n");
}

