#ifndef __EXPOTREEEXT_H__
#define __EXPOTREEEXT_H__

#include <math.h>
#include <R.h>
#include <Rdefines.h>

void expotree_ext(double* K, double* beta, double* mu,
                  double* psi, int* parVecLen, int* ki, int* n, 
                  double* times, int* ttypes, double* p, double* t0, 
                  int* info, int* estNorm, int* vflag, int* rescale);

SEXP expoTreeEvalExt(SEXP K, SEXP lambda, SEXP mu, SEXP psi, SEXP rho,
                     SEXP times, SEXP ttypes, SEXP add_par);

#endif
