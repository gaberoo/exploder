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

#ifndef __simulate_trees_h__
#define __simulate_trees_h__

#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <cmath>
#include <typeinfo>
#include <string>
#include <sstream>
#include <algorithm>
#include <getopt.h>
using namespace std;

#ifdef FAKER
#include <GSLRng.h>
#else
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#endif

typedef enum { BD = 0, SIS = 1, SIR = 2 } SimType; 

class Pars {
public:
  Pars() : N(-1), K(-1.0), beta(1.0), mu(0.1), psi(0.1), t(0.0), lfun(SIS) {}
  virtual ~Pars() {}

  double& b() { return beta; }
  double& u() { return mu; }
  double& s() { return psi; }

  int N;
  double K;
  double beta;
  double mu;
  double psi;
  double t;

  SimType lfun;
};

// ===========================================================================

inline double lambda(int I, const Pars& pars) { 
  return fmax(0.0,pars.beta*(1.-I/pars.K));
}

inline double lambdaSIR(int I, int S, const Pars& pars) { 
  return fmax(0.0,pars.beta*S/pars.K);
}

// ===========================================================================

struct Individual {
  Individual(double it = 0.0, int p = 0) 
    : itime(it), dtime(-1.0), stime(-1.0), parent(p), id(-1)
  {}

  virtual ~Individual() {}

  double itime;  // infection time
  double dtime;  // death time
  double stime;  // sample time
  int parent;    // parent individual
  int id;        // my id
  vector<int> children;
};

// ===========================================================================

int getBranchingTimes(vector<double>& times, int subclade, 
    const vector<Individual>& pop);

string newickString(int subclade, const vector<Individual>& pop, 
    double t, bool sampled = false);

double sim_trees(const Pars& pars, vector<Individual>& pop, list<int>& inf, vector<int>& samples,
    int max_samples, double max_time, int min_outbreak, int resume = 0, double curTime = 0.0);

void to_array(double t, vector<Individual>& pop, vector<int>& samples, 
    size_t len, double* x, int* xtype);

void inf_times(double t, vector<Individual>& pop, 
    size_t len, double* x0, double* x1, int* xtype, 
    int* id, int* parent);

#endif // __simulate_trees_h__

