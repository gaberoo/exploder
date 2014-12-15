#include "expoTree.h"

SEXP expoTreeEval(SEXP parameters, SEXP times, SEXP ttypes, SEXP add_par) {
  PROTECT(parameters = AS_NUMERIC(parameters));
  PROTECT(times = AS_NUMERIC(times));
  PROTECT(ttypes = AS_INTEGER(ttypes));
  PROTECT(add_par = AS_INTEGER(add_par));

  Rboolean parMat = isMatrix(parameters);
  int parVecLen = 1;
  int parVecCol = 1;
  int quit = 0;
  int surv = INTEGER(add_par)[0];
  double zeroTol = 1e-12;

  int vf = (LENGTH(add_par) > 1) ? INTEGER(add_par)[1] : 0;
  if (vf) Rprintf("Verbose mode is on.\n");

  int rs = (LENGTH(add_par) > 2) ? (INTEGER(add_par)[2]>0) : 1;
  if (rs && vf) Rprintf("Rescaling is on.\n");

  int extant = (LENGTH(add_par) > 3) ? INTEGER(add_par)[3] : 0;
  if (vf) Rprintf("%d lineages extant at the root.\n",extant);

  int root = (LENGTH(add_par) > 4) ? INTEGER(add_par)[4] : 0;
  if (root && vf) Rprintf("Root correction requested.\n");

  int est_norm = (LENGTH(add_par) > 5) ? INTEGER(add_par)[5] : 0;
  if (est_norm && vf) Rprintf("Forcing estimation of matrix norm.\n");

  // check dimensions
  if (parMat) {
    if (vf) Rprintf("Checking matrix dimensions: ");
    SEXP parDim = GET_DIM(parameters);
    parVecLen = INTEGER(parDim)[0];
    parVecCol = INTEGER(parDim)[1];
    if (vf) Rprintf("nrow = %d, ncol = %d\n",parVecLen,parVecCol);
    if (parVecCol < 5 || parVecLen < 1) quit = 1;
  } else {
    if (LENGTH(parameters) < 5) {
      if (vf) Rprintf("Error: parameter vector must have a minimum length of 5!\n");
      quit = 1;
    }
  }

  if (quit) {
    // not enough columns
    if (vf) Rprintf("Exiting prematurely.\n");
    SEXP p;
    PROTECT(p = NEW_NUMERIC(1));
    REAL(p)[0] = R_NegInf;
    UNPROTECT(5);
    return p;
  }

  int i;

  double* pars = NUMERIC_POINTER(parameters);

  double* N    = pars;
  double* beta = pars+parVecLen;
  double* mu   = pars+2*parVecLen;
  double* psi  = pars+3*parVecLen;
  double* drho = pars+4*parVecLen;

  /* only use the first rho entry */
  double rho = drho[0];

  /* convert N reals to integers */
  int Nmax = 0;
  for (i = 0; i < parVecLen; ++i) {
    if (Nmax < ceil(N[i])) Nmax = (int) ceil(N[i]);
  }
  if (vf) Rprintf("Maximum dimension: N = %d\n",Nmax);

  int SI = 1;

  double* ptimes  = NUMERIC_POINTER(times);
  int*    pttypes = INTEGER_POINTER(ttypes);
  int     nt      = LENGTH(times);

  int maxExtant = 0;

  int ki = 0;

  SEXP p;
  PROTECT(p = NEW_NUMERIC(Nmax+1));  /* +1 for p(0) */
  double* p0 = NUMERIC_POINTER(p);

  double t0 = 0.0;
  double scale = 0.0;

  for (i = nt-1; i >= 0; --i) {
    switch (pttypes[i]) {
      case 1:
      case 11:
        ++extant;
        break;
      case 0:
      case 2:
      case 4:
      case 10:
        --extant;
        break;
      case 3:
      default:
        break;
    }
    if (maxExtant < extant) maxExtant = extant;
  }
  if (vf) {
    Rprintf("%d extant lineages at the present; %d maximum.\n",
            extant,maxExtant);
  }

  int goodParams = 1;
  for (i = 0; i < parVecLen && goodParams; ++i) {
    if (beta[i] <= 0.0 || mu[i] < 0.0 || psi[i] < 0.0) goodParams = 0;
  }

  if (! goodParams || rho < 0.0 || rho > 1.0) {
    if (vf) {
      Rprintf("Illegal parameters:\n");
      for (i = 0; i < parVecLen; ++i) {
        Rprintf("%f %f %f %f %f\n",N[i],beta[i],mu[i],psi[i],drho[i]);
      }
    }
    for (i = 0; i <= Nmax; ++i) p0[i] = R_NegInf;
  } else {
    if (surv) {
      p0[0] = 1.0;
      for (i = 1; i <= Nmax; ++i) {
        p0[i] = (rho <= 0.0 || rho >= 1.0) ? 0.0 : pow(1.-rho,i);
      }
      rExpoTree(N,&ki,beta,mu,psi,&nt,&parVecLen,ptimes,pttypes,
                p0,&t0,&SI,&est_norm,&vf,&rs);
      for (i = 0; i <= N[parVecLen-1]; ++i) {
        if (p0[i] >= zeroTol) {
          Rprintf("Log survival probability is non-negative! p(%d) = %g\n",i,p0[i]);
        }
        p0[i] = (p0[i] < 0.0) ? log(1.-exp(p0[i])) : R_NegInf;
      }
    } else {
      // set initial value of p
      if (extant == 0) {
        p0[0] = 0.0;
        for (i = 1; i <= Nmax; ++i) p0[i] = psi[0];
        ki = 1;
        t0 = ptimes[0];
        ptimes = ptimes+1;
        pttypes = pttypes+1;
        --nt;
      } else {
        ki = extant;
        p0[0] = 0.0;
        scale = extant*log(rho);
        for (i = 1; i <= Nmax; ++i) {
          if (i < extant) p0[i] = 0.0;
          else p0[i] = pow(1.-rho,i-extant);
        }
      }
      rExpoTree(N,&ki,beta,mu,psi,&nt,&parVecLen,ptimes,pttypes,
                p0,&t0,&SI,&est_norm,&vf,&rs);
      for (i = 0; i <= N[parVecLen-1]; ++i) {
        p0[i] += scale;
        if (root) {
          p0[i] -= M_LN2 + log(beta[parVecLen-1]) 
                   + log(1.-1./N[parVecLen-1]);
        }
      }
    }
  }

  UNPROTECT(5);
  return p;
}

/****************************************************************************/

double bdss_q(double t, double c1, double c2) {
  double q = (1-c2)*exp(-t*c1/2.0) + ((1+c2)*exp(t*c1/2));
  return 0.25*q*q;
}

/****************************************************************************/

SEXP infTreeEval(SEXP parameters, SEXP times, SEXP ttypes, SEXP add_par) {
  PROTECT(parameters = AS_NUMERIC(parameters));
  PROTECT(times = AS_NUMERIC(times));
  PROTECT(ttypes = AS_INTEGER(ttypes));
  PROTECT(add_par = AS_INTEGER(add_par));

  double* pars = NUMERIC_POINTER(parameters);
  double beta = pars[0];
  double mu   = pars[1];
  double psi  = pars[2];
  double rho  = pars[3];

  double* ptimes = NUMERIC_POINTER(times);
  int* pttypes = INTEGER_POINTER(ttypes);
  int nt = LENGTH(times);
  int maxExtant = 0;

  int surv = INTEGER(add_par)[0];

  int vf = (LENGTH(add_par) > 1) ? INTEGER(add_par)[1] : 0;
  if (vf) Rprintf("Verbose mode is on.\n");

  int extant = (LENGTH(add_par) > 2) ? INTEGER(add_par)[2] : 0;
  if (vf) Rprintf("%d lineages extant at the root.\n",extant);

  SEXP fx;
  PROTECT(fx = NEW_NUMERIC(1));
  double* lik = NUMERIC_POINTER(fx);

  int root = (ptimes[nt-1] == ptimes[nt-2]) ? 1 : 0;

  double t0 = 0.0;
  int i;

  for (i = nt-1; i >= 0; --i) {
    switch (pttypes[i]) {
      case 1:
        ++extant;
        break;
      case 0:
      case 2:
      case 4:
        --extant;
        break;
      default:
        break;
    }
    if (maxExtant < extant) maxExtant = extant;
  }

  if (beta <= 0.0 || mu < 0.0 || psi < 0.0 || rho < 0.0 || rho > 1.0) {
    if (vf) Rprintf("Some parameters are out of range.");
    *lik = R_NegInf;
  } else {
    double c1;
    double c2;
    double p0;
    double mp = mu+psi;

    c1 = beta-mp;
    c1 = sqrt(c1*c1 + 4*beta*psi);
    c2 = -(beta-mp-2*beta*rho)/c1;
    *lik = -log(2*beta);

    if (surv) {
      p0 = exp(-c1*ptimes[nt-1])*(1.-c2);
      p0 = beta+mp+c1*(p0-(1.+c2))/(p0+1.+c2);
      p0 = p0/(2.*beta);
      *lik -= log(1.-p0);
    }

    if (extant > 0) *lik += extant*log(rho);
    for (int i = 0; i < nt; ++i) {
      switch (pttypes[i]) {
        case 1:
          *lik += log(2*beta/bdss_q(ptimes[i],c1,c2));
          break;

        case 0:
          *lik += log(psi*bdss_q(ptimes[i],c1,c2));
          break;

        default:
          break;
      }
    }
  }

  UNPROTECT(5);
  return fx;
}


