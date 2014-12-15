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

#include <gsl/gsl_rng.h>
gsl_rng* globalRng;

#define FAKER
#include "R.h"

#include "shared/simulate_trees.h"
#include <gsl/gsl_math.h>

/****************************************************************************/

void printHelp() {
  printf("sim_trees : Simulate transmission trees under epidemiological models.\n\n");
  printf("usage: sim_trees \n");
  printf("-N <population size> ");
  printf("-b <beta> -m <mu> -p <psi> [-tnAh] [-s <samples>] [-O <infecteds>] [-S <seed>]\n\n");
  printf("  N : total population size\n");
  printf("  b : infection rate 'beta'\n");
  printf("  u : recovery rate 'mu'\n");
  printf("  s : sampling rate 'psi'\n");
  printf("  t : model type\n");
  printf("       0 = constant rate\n");
  printf("       1 = SIS\n");
  printf("       2 = SIR\n");
  printf("       3 = SEIR\n");
  printf("  d : output mode\n");
  printf("       1 = full tree in Nexus format\n");
  printf("       2 = branching times (default)\n");
  printf("       3 = epidemic curve\n");
  printf("       0 = all\n");
  printf("  n : output tree in Newick format\n");
  printf("  A : output both the branching times and the tree\n");
  printf("  S : set the random number seed\n");
  printf("  x : exit simulation when s individuals have been sampled\n");
  printf("  O : minimum number of infecteds considered an outbreak\n");
  printf("  h : output this help message\n");
}

/****************************************************************************/

class Event {
  public:
    Event(int a = 0, double b = 0.0) : type(a), time(b) {}
    virtual ~Event() {}
    int type;
    double time;
};

int compare_event(const void * a, const void * b) {
  Event* aa = (Event*) a;
  Event* bb = (Event*) b;
  if (aa->time < bb->time) return -1;
  else if (aa->time > bb->time) return 1;
  else return 0;
}

int main(int argc, char* argv[]) {
  globalRng = gsl_rng_alloc(gsl_rng_taus2);

  Pars pars;
  pars.N    = -1;
  pars.beta = -1.0;
  pars.mu   = -1e-3;
  pars.psi  = -1e-3;
  pars.lfun = SIS;

  int c;
  opterr = 0;
  int output = 2;
  int seed = time(NULL);
  int max_samples = 100;
  int min_outbreak = 20;
  double max_time = INFINITY;

  while ((c = getopt(argc,argv,"N:b:u:s:t:d:S:x:O:hT:")) != -1) {
    switch (c) {
      case 'N': 
        pars.N = atoi(optarg); 
        pars.K = 1.0*pars.N;
        break;
      case 'b': pars.beta = atof(optarg); break;
      case 'u': pars.mu = atof(optarg); break;
      case 's': pars.psi = atof(optarg); break;
      case 't': pars.lfun = (SimType) atoi(optarg); break;
      case 'd': output = atoi(optarg); break;
      case 'S': seed += atoi(optarg); break;
      case 'x': max_samples = atoi(optarg); break;
      case 'O': min_outbreak = atoi(optarg); break;
      case 'T': max_time = atof(optarg); break;
      case 'h': printHelp(); return 0;
    }
  }

  min_outbreak = GSL_MIN(max_samples,min_outbreak);

  if (pars.beta < 0.0 || pars.mu < 0.0 || pars.psi < 0.0) {
    fprintf(stderr,"Please specify beta, mu and psi.\n\n");
    return 0;
  }

  gsl_rng_set(globalRng,seed);

  vector<Individual> pop;
  list<int> inf;
  vector<int> samples;

  if (max_time < 0) max_time = R_PosInf;

  double t = sim_trees(pars,pop,inf,samples,max_samples,
                       max_time,min_outbreak);

  if (t == 0.0) cerr << "Couldn't simulate outbreak." << endl;

  if (output == 1 || output == 0) {
    pars.t = t;
    cout << "#NEXUS" << endl;
    cout << "begin trees;" << endl;
    cout << "tree 'simTree' = "
         << newickString(0,pop,0.0,true) << ";" << endl;
    cout << "end;" << endl;
  }

  if (output == 2 || output == 0) {
    cout << "#> --- TIMES ---" << endl;
    printf("#  N = %d\n",pars.N);
    printf("#  b = %g\n",pars.beta);
    printf("#  m = %g\n",pars.mu);
    printf("#  p = %g\n",pars.psi);
    printf("#  S = %d\n",seed);
    printf("#  s = %d\n",max_samples);
    printf("#  O = %d\n",min_outbreak);

    int maxlen = max_samples*2;
    double* x = new double[maxlen];
    int* xt = new int[maxlen];
    to_array(t,pop,samples,maxlen,x,xt);

    vector<double> times;
    getBranchingTimes(times,0,pop);
    times.push_back(0.0);

    for (int i(0); i < times.size(); ++i) {
      cout << t-times[i] << " 1 0" << endl;
    }

    for (int i(0); i < samples.size(); ++i) {
      cout << t-pop[samples[i]].stime << " 0 " << samples[i] << endl;
    }

    delete[] x;
    delete[] xt;
  } 

  if (output == 3 || output == 0) {
    cout << "#> --- CUMULATIVE INFECTIONS ---" << endl;
    int maxlen = max_samples*2;
    int itlen = pop.size();

    vector<double> times(maxlen,0.0);
    vector<int>    ttypes(maxlen,0);
    vector<double> itimes(itlen,0.0);
    vector<double> dtimes(itlen,0.0);
    vector<int>    dtypes(itlen,0);
    vector<int>    id(itlen,0);
    vector<int>    parent(itlen,0);

    to_array(t,pop,samples,maxlen,times.data(),ttypes.data());
    inf_times(t,pop,itlen,itimes.data(),dtimes.data(),
              dtypes.data(),id.data(),parent.data());

    vector<Event> events;
    for (int i(0); i < itlen; ++i) {
      events.push_back(Event(1,t-itimes[i]));
      events.push_back(Event(2,t-dtimes[i]));
    }

    qsort(events.data(),events.size(),sizeof(Event),compare_event);

    int cumsum = 0;
    for (int i(0); i < events.size(); ++i) {
      switch (events[i].type) {
        case 1: ++cumsum; break;
        case 2: --cumsum; break;
        default: break;
      }
      cout << events[i].time << " " << cumsum << endl;
    }
  }

  gsl_rng_free(globalRng);

  return 0;
}

/****************************************************************************/


