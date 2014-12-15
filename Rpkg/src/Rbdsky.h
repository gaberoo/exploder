#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

double bdsky_eval(const double* beta, const double* mu, const double* psi, const double* rho,
                  int nshifts, const double* shifts, int n, const double* times, const int* ttypes, 
                  int survival, int vflag);

SEXP bdSkyEval(SEXP parameters, SEXP times, SEXP ttypes, SEXP shifts, SEXP add_par);

