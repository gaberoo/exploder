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

#ifndef __EXPOTREE_H__
#define __EXPOTREE_H__

#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_permute_int.h>

extern "C" {
  #include "shared/rexpotree.h"

  void sir_expotree(double* RN, int* Rki, double* Rbeta, double* Rmu,
      double* Rpsi, int* Rn, int* parVecLen, 
      double* times, int* ttypes, double* p, double* t0, 
      int* info, int* estimateNorm, int* Rvflag, int* Rrescale);
}

#ifdef OPENCL
extern "C" {
  void clrExpoTree(double* RN, int* Rki, double* Rbeta, double* Rmu,
      double* Rpsi, int* Rn, int* parVecLen, double* times, 
      int* ttypes, double* p, double* t0, int* info, 
      int* estimateNorm, int* Rvflag, int* Rrescale);
}
#endif

void readTimes(string fn, vector<double>& times, vector<int>& ttypes,
               int& extant, int& maxExtant);

/*
 * N         : maximum population size
 * beta      : infection/speciation rate
 * mu        : recovery/extinction rate
 * psi       : sampling rate
 * rho       : initial sampling probability
 * times     : vector of event times
 * ttypes    : vector of event types
 * extant    : number of extant species at the present
 * SI        : model type (SI = 1: density-dependence, SI = 0: no
 *             density-dependence)
 * vf        : verbosity level
 * rs        : rescale probability vector
 */

/****************************************************************************/

double expoTreeEval(int parVecLen, double* K, 
    double* beta, double* mu, double* psi, double rho, 
    vector<double>& times, vector<int>& ttypes, 
    int extant, int est_norm, int vflag, int rescale,
    int nroot = 0);

inline double expoTreeEval(vector<double>& K, vector<double>& beta, 
    vector<double>& mu, vector<double>& psi, double rho, 
    vector<double>& times, vector<int>& ttypes, 
    int extant, int estimateNorm, int vflag, int rescale,
    int nroot = 0) 
{
  return expoTreeEval((int) K.size(),K.data(),beta.data(),mu.data(),
      psi.data(),rho,times,ttypes,extant,estimateNorm,vflag,rescale,
      nroot);
}

/****************************************************************************/

double expoTreeSurvival(vector<double>& K, vector<double>& beta, 
    vector<double>& mu, vector<double>& psi, double rho,
    vector<double>& times, vector<int>& ttypes, 
    int extant, int SI, int vf, int rs);

// wrapper functions (for backwards compatibility)

inline double expoTreeEval(double K, double beta, double mu, double psi, 
    double rho, vector<double>& times, vector<int>& ttypes, 
    int extant, int SI, int vf, int rs) {
  vector<double> vK(1,K);
  vector<double> vBeta(1,beta);
  vector<double> vMu(1,mu);
  vector<double> vPsi(1,psi);
  double lik = expoTreeEval(vK,vBeta,vMu,vPsi,rho,times,ttypes,extant,SI,vf,rs);
  return lik;
}

inline double expoTreeSurvival(double K, double beta, double mu, double psi, 
    double rho, vector<double>& times, vector<int>& ttypes, 
    int extant, int SI, int vf, int rs) {
  vector<double> vK(1,K);
  vector<double> vBeta(1,beta);
  vector<double> vMu(1,mu);
  vector<double> vPsi(1,psi);
  return expoTreeSurvival(vK,vBeta,vMu,vPsi,rho,times,ttypes,extant,SI,vf,rs);
}

inline double expoTreeSurvival(double K, double beta, double mu, double psi, 
    double rho, double torig, int extant, int SI, int vf, int rs) {
  vector<double> vK(1,K);
  vector<double> vBeta(1,beta);
  vector<double> vMu(1,mu);
  vector<double> vPsi(1,psi);
  vector<double> times(1,torig);
  vector<int> ttypes(1,1);
  return expoTreeSurvival(vK,vBeta,vMu,vPsi,rho,times,ttypes,extant,SI,vf,rs);
}


// ===========================================================================

inline double bdss_q(double t, double c1, double c2) {
  double q((1-c2)*exp(-t*c1/2.0) + ((1+c2)*exp(t*c1/2)));
  return 0.25*q*q;
}

double infTreeSurvival(double beta, double mu, double psi, double rho, 
    double torig);

double infTreeEval(double beta, double mu, double psi, double rho,
    const vector<double>& times, const vector<int>& ttypes, 
    int extant, int SI, int survival, int vf);

int check_params(int n, double* N, double* beta, double* mu, 
                 double* psi, double rho);

#endif
