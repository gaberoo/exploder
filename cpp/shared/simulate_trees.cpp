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

#include "simulate_trees.h"

// ===========================================================================

int getBranchingTimes(vector<double>& times, int subclade, 
    const vector<Individual>& pop) 
{
  int ret = 0;
  // check if the leaf is sampled
  if (pop[subclade].stime >= 0.0) ret = 1;
  // check if at least one of the children is sampled
  for (int i(pop[subclade].children.size()-1); i >= 0; --i) {
    if (getBranchingTimes(times,pop[subclade].children[i],pop)) {
      if (ret) {
        times.push_back(pop[pop[subclade].children[i]].itime);
      } else {
        ret = 1;
      }
    }
  }
  return ret;
}

// ===========================================================================

string newickString(int subclade, const vector<Individual>& pop, 
    double t, bool sampled) 
{
  ostringstream str;
  double dtime(0.0);
  if (pop[subclade].stime >= 0.0) dtime = pop[subclade].stime;
  else if (pop[subclade].dtime >= 0.0) dtime = pop[subclade].dtime;
  else dtime = t;
  if (pop[subclade].children.size() == 0) {
    str << "'" << subclade;
    if (pop[subclade].stime >= 0.0) str << "+";
    str << "':";
    str << dtime-pop[subclade].itime;
  } else {
    t = dtime;
    for (size_t i(0); i < pop[subclade].children.size(); ++i) str << "(";
    str << "'" << subclade;
    if (pop[subclade].stime >= 0.0) str << "+";
    str << "'";
    for (int i(pop[subclade].children.size()-1); i >= 0; --i) {
      str << ":" << t-pop[pop[subclade].children[i]].itime;
      str << "," << newickString(pop[subclade].children[i],pop,t,sampled) << ")";
      // str << "'" << subclade << "->" << pop[subclade].children[i] << "'";
      t = pop[pop[subclade].children[i]].itime;
    }
    str << ":" << t-pop[subclade].itime;
  }
  return str.str();
}

// ===========================================================================

double sim_trees(const Pars& pars, 
                 vector<Individual>& pop,   /* whole population */
                 list<int>& inf,            /* pointers to infecteds */
                 vector<int>& samples,      /* pointers to sampleds */
                 int max_samples,           /* stop after this # samples */
                 double max_time,           /* stop after this time */
                 int min_outbreak,          /* minimum infections that counts as an outbreak */
                 int resume,                /* (flag) resume previous run */
                 double curTime)            /* current time */
{
  /* Reset RNG state (only for R mode) */
  GetRNGstate();

  double t(0.0);              /* running time */
  double totalRate(0.0);      /* sum of all rates */
  double nextEventTime(0.0);  /* next event happening at time */

  int S(pars.N-1); /* susceptibles */
  int I(1);        /* infecteds */
  int J(0);        /* sampled recovereds */

  int node(0);

  list<int>::iterator infit;

  int max_tries = 1000;
  int tries = 0;

  vector<Individual> popBackup;
  if (resume) popBackup = pop;

  while (tries++ < max_tries) { /* retry is early extinction */
    /* clear all lists... */
    inf.clear();
    samples.clear();

    /* ...and reset times and rates */
    t = curTime;
    nextEventTime = 0.0;

    /* resume a previous run (TODO for SIR) */
    if (resume) {
      pop = popBackup;
      I = 0;
      J = 0;
      for (size_t i(0); i < pop.size(); ++i) {
        if (pop[i].dtime < 0.0 && pop[i].stime < 0.0) {
          inf.push_back(i);
          ++I;
        } else if (pop[i].stime >= 0.0) {
          samples.push_back(i);
          ++J;
        }
      }
    } else {
      pop.clear();
      // infect initial
      pop.push_back(Individual(0.0,-1));
      pop.back().id = 0;
      inf.push_back(0);
      I = 1;
      J = 0;
    }

    /* run while 
     *   (i)   there are still infecteds left; 
     *   (ii)  the number of sampled individuals is not yet reached;
     *   (iii) the maximum time is not yet reached */
    while ((I+S) > 0 
           && samples.size() < (size_t) max_samples 
           && t < max_time) 
    {
      /* if psi = 0 (no samples taken), then the number of infecteds is the
       * threshold for the number of samples */
      if (pars.psi <= 0.0 && I >= max_samples) {
        /* in this case, add all infecteds to the sampled list */
        for (infit = inf.begin(); infit != inf.end(); ++infit) {
          pop[*infit].stime = t;
          samples.push_back(*infit);
        }
        break;
      }

      double l = 0.0;
      /* calculate lambda */
      switch (pars.lfun) {
        case SIS: l = lambda(I,pars); break;
        case SIR: l = lambdaSIR(I,S,pars); break;
        case BD:  
        default: l = pars.beta; break;
      }

      /* sum up to total rate */
      totalRate = l*I + (pars.mu+pars.psi)*I;

      if (totalRate <= 0.0) break;

      /* sample next event time */
      nextEventTime = t-log(unif_rand())/totalRate;

      /* pick random individual to whom the event happens */
      node = floor(unif_rand()*I);
      infit = inf.begin();
      advance(infit,node);

      /* pick event */
      double r = unif_rand();

      if (r < l*I/totalRate) {
        /* branching/infection event */
        pop.push_back(Individual());
        pop.back().itime = nextEventTime;
        pop.back().parent = *infit;
        pop.back().id = pop.size()-1;
        inf.push_back(pop.size()-1);
        pop[*infit].children.push_back(pop.size()-1);
        ++I;
        --S;
      } else {
        /* removal/recovery event */
        if (r < (l*I+pars.mu*I)/totalRate) {
          pop[*infit].dtime = nextEventTime;
        } else {
          pop[*infit].stime = nextEventTime;
          samples.push_back(*infit);
          ++J;
        }
        inf.erase(infit);
        --I;
        if (pars.lfun == SIS) ++S;
      }

      t = nextEventTime;
    }

    if (pars.lfun == SIR || (int) samples.size() >= min_outbreak) break;
  }

  PutRNGstate();

  return t;
}

// ===========================================================================

void to_array(double t, vector<Individual>& pop, vector<int>& samples, 
    size_t len, double* x, int* xtype) {
  vector<double> times;
  getBranchingTimes(times,0,pop);
  times.push_back(0.0);
  int j(0);
  if (len >= times.size() + samples.size()) {
    for (size_t i(0); i < times.size(); ++i) {
      x[j] = t-times[i];
      xtype[j] = 1;
      ++j;
    }
    for (size_t i(0); i < samples.size(); ++i) {
      x[j] = t-pop[samples[i]].stime;
      xtype[j] = 0;
      ++j;
    }
  }
}

// ===========================================================================

void inf_times(double t, vector<Individual>& pop, 
    size_t len, double* x0, double* x1, int* xtype, 
    int* id, int* parent) 
{
  // loop through all individuals of the population
  // while (ind != pop.end() && i < len) {
  for (size_t i(0); i < pop.size(); ++i) {
    // Rprintf("[%4d] %g\n",pop[i].id,pop[i].itime);
    // if (pop[i].itime >= 0) {
      x0[i] = t-pop[i].itime;
      id[i] = pop[i].id;
      parent[i] = pop[i].parent;
      if (pop[i].stime >= 0) {
        x1[i] = t-pop[i].stime;
        xtype[i] = 0;
      } else if (pop[i].dtime >= 0) {
        x1[i] = t-pop[i].dtime;
        xtype[i] = -1;
      } else {
        x1[i] = 0.0;
        xtype[i] = 1;
      }
    // }
  }
}

// ===========================================================================

