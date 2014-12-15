#ifndef __EXPOTREESIR_H__
#define __EXPOTREESIR_H__

#include "expoTree.h"

extern "C" {
  int sir_index_N(int I, int R, int N);

  void sir_expotree(double* RN, int* Rki, double* Rbeta, double* Rmu,
      double* Rpsi, int* Rn, int* parVecLen, 
      double* times, int* ttypes, double* p, double* t0, 
      int* info, int* estimateNorm, int* Rvflag, int* Rrescale);
}

double expoTreeSIR(int parVecLen, double* K, 
                    double* beta, double* mu, double* psi, 
                    double rho, vector<double>& times, vector<int>& ttypes, 
                    int extant, int est_norm, int vflag, int rescale,
                    int nroot);

double expoTreeSIRSurvival(int parVecLen, double* K, double* beta, double* mu, 
                           double* psi, double rho, vector<double>& times, 
                           vector<int>& ttypes, int extant, int est_norm, 
                           int vflag, int rescale, int nroot);

#endif
