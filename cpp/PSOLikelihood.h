#ifndef _PSOLIK_H_
#define _PSOLIK_H_

#include "pso/Point.h"
#include "expoTree.h"
#include "Forest.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

// ===========================================================================

typedef struct {
  int extant;
  int maxExtant;
  int SImodel;
  int vflag;
  int rescale;
  double lockedRatio;
  Forest* trees;
} PSOParams;

// ===========================================================================

double likelihood(const PSO::Point& x, void* params) {
  PSOParams par = *(PSOParams*) params;

  double fx;
  int    N    = (int) round(x[4]); 
  double beta = x[0];
  double mu   = (par.lockedRatio >= 0.0) ? x[2]*(1./par.lockedRatio-1.) : x[1];
  double psi  = x[2];
  double rho  = x[3];
  int vflag = (par.vflag > 0) ? par.vflag-1 : vflag;

  vector<double> vars(x);
  vars[1] = mu;

  if (N <= par.maxExtant || beta <= 0.0 || mu < 0.0 || psi < 0.0 
      || rho < 0.0 || rho > 1.0) {
    fx = -INFINITY;
  } else {
    if (par.SImodel)
      fx = par.trees->jointLikelihood(vars,par.SImodel,vflag,par.rescale);
    else
      fx = par.trees->jointInfLik(vars,par.SImodel,par.vflag,par.rescale);
  }

  return fx;
}

// ===========================================================================

double likelihood_inf(const PSO::Point& x, void* params) {
  PSOParams par = *(PSOParams*) params;

  double fx = 0.0;

  double beta = x[0];
  double mu   = (par.lockedRatio >= 0.0) ? x[2]*(1./par.lockedRatio-1.) : x[1];
  double psi  = x[2];
  double rho  = x[3];

  vector<double> vars(x);
  vars[1] = mu;

  if (beta <= 0.0 || mu < 0.0 || psi < 0.0 || rho < 0.0 || rho > 1.0) {
    fx = -INFINITY;
  } else {
    fx = par.trees->jointInfLik(vars,par.SImodel,par.vflag,par.rescale);
  }

  return fx;
}

#endif
