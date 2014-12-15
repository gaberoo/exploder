#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "ExpoTreePars.h"

typedef struct _pars {
  int lfun;
  double N;
  double beta;
  double mu;
} Pars;

int SIR (double t, const double y[], double f[], void* pp) {
  Pars p = *(Pars*) pp;
  f[0] = p.beta*(p.N-y[0]-y[1])*y[0]/p.N - p.mu*y[0];
  f[1] = p.mu*y[0];
  return GSL_SUCCESS;
}

int SIR_jac (double t, const double y[], double* dfdy, double dfdt[], void* pp) {
  Pars p = *(Pars*) pp;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy,2,2);
  gsl_matrix* m = &dfdy_mat.matrix; 
  gsl_matrix_set(m, 0, 0, -p.mu + p.beta*(p.N-y[1]-2*y[0])/p.N);
  gsl_matrix_set(m, 0, 1, -p.beta/p.N);
  gsl_matrix_set(m, 1, 0, p.mu);
  gsl_matrix_set(m, 1, 1, 0.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

int main(int argc, char* argv[]) {
  int nroot = 1;
  int vflag;
  int model_type = 1;

  double maxTime = 0.0;

  vector<double> K;
  vector<double> beta;
  vector<double> mu;
  vector<double> psi;
  vector<double> rShift;

  int c;
  opterr = 0;
  while ((c = getopt(argc,argv,"n:N:vb:u:s:S:t:T:")) != -1) {
    switch (c) {
      case 'n': nroot = atoi(optarg); break;
      case 'N': K.push_back(atof(optarg)); break;
      case 'v': ++vflag; break;
      case 'b': beta.push_back(atof(optarg)); break;
      case 'u': mu.push_back(atof(optarg)); break;
      case 's': psi.push_back(atof(optarg)); break;
      case 'S': rShift.push_back(atof(optarg)); break;
      case 't': model_type = atoi(optarg); break;
      case 'T': maxTime = atof(optarg); break;
      default: break;
    }
  }

  // =========================================================================

  Pars p = { 2, 10.0, 1.0, 0.1 };

  gsl_odeiv2_system sys;

  switch (model_type) {
    case 2:
      sys.function = SIR;
      sys.jacobian = SIR_jac;
      sys.dimension = 2;
      sys.params = &p;
      break;

    case 1:
    default:
      break;
  }

  gsl_odeiv2_driver * d = 
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45,
                                   1e-6, 1e-6, 0.0);
  double y[2] = { 1.0*nroot, 0.0 };

  int nshifts = rShift.size();

  double curTime = 0.0;
  double nextTime = 0.0;
  double ti;
  double dt = 1.0;
  int status;
  int j;
  for (j = nshifts; j >= 0; --j) {
    nextTime = (j>0) ? maxTime - rShift[j-1] : maxTime;
    p.N = K[j];
    p.beta = beta[j];
    p.mu = (mu[j]+psi[j]);
    while (curTime < nextTime) {
      ti = GSL_MIN(curTime+dt,nextTime);
      status = gsl_odeiv2_driver_apply (d, &curTime, ti, y);
      if (status != GSL_SUCCESS) {
        printf("error, return value=%d\n", status);
        exit(1);
      }
      printf ("%.5e %.5e %.5e\n", curTime, y[0], y[1]);
    }
  }

  gsl_odeiv2_driver_free (d);
  return 0;
}
