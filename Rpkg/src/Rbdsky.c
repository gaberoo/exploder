#include "Rbdsky.h"

double bdsky_eval(const double* beta, const double* mu, const double* psi, const double* rho,
                  int nshifts, const double* shifts, int n, const double* times, const int* ttypes, 
                  int survival, int vflag);

SEXP bdSkyEval(SEXP parameters, SEXP times, SEXP ttypes, SEXP shifts, SEXP add_par) {
  PROTECT(parameters = AS_NUMERIC(parameters));
  PROTECT(times = AS_NUMERIC(times));
  PROTECT(ttypes = AS_INTEGER(ttypes));
  PROTECT(shifts = AS_NUMERIC(shifts));
  PROTECT(add_par = AS_INTEGER(add_par));

  Rboolean parMat = isMatrix(parameters);
  int parVecLen = 1;
  int parVecCol = 1;
  int quit = 0;
  int surv = INTEGER(add_par)[0];
  double zeroTol = 1e-12;

  int vf = (LENGTH(add_par) > 1) ? INTEGER(add_par)[1] : 0;
  if (vf) Rprintf("Verbose mode is on.\n");

  int nroot = (LENGTH(add_par) > 2) ? INTEGER(add_par)[2] : 0;
  if (vf) Rprintf("Root lineages = %d.\n",nroot);

  // check dimensions
  if (parMat) {
    if (vf) Rprintf("Checking matrix dimensions: ");
    SEXP parDim = GET_DIM(parameters);
    parVecLen = INTEGER(parDim)[0];
    parVecCol = INTEGER(parDim)[1];
    if (vf) Rprintf("nrow = %d, ncol = %d\n",parVecLen,parVecCol);
    if (parVecCol < 4 || parVecLen < 1) quit = 1;
  } else {
    if (LENGTH(parameters) < 4) {
      if (vf) Rprintf("Error: parameter vector must have a minimum length of 4!\n");
      quit = 1;
    }
  }

  if (quit) {
    // not enough columns
    if (vf) Rprintf("Exiting prematurely.\n");
    SEXP p;
    PROTECT(p = NEW_NUMERIC(1));
    REAL(p)[0] = R_NegInf;
    UNPROTECT(6);
    return p;
  }

  int i;

  double* pars = NUMERIC_POINTER(parameters);

  double* beta = pars;
  double* mu   = pars+  parVecLen;
  double* psi  = pars+2*parVecLen;
  double* rho  = pars+3*parVecLen;

  double* ptimes  = NUMERIC_POINTER(times);
  int*    pttypes = INTEGER_POINTER(ttypes);
  int     nt      = LENGTH(times);

  double* pshifts = NUMERIC_POINTER(shifts);
  int     nshifts = LENGTH(shifts);

  SEXP p;
  PROTECT(p = NEW_NUMERIC(1));

  REAL(p)[0] = bdsky_eval(beta,mu,psi,rho,nshifts,pshifts,nt,ptimes,pttypes,surv,vf);

  UNPROTECT(6);
  return p;
}

