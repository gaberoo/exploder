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

#include "expoTreeExt.h"

SEXP expoTreeEvalExt(SEXP K, SEXP lambda, SEXP mu, SEXP psi, SEXP rho,
                     SEXP times, SEXP ttypes, SEXP add_par) 
{
  PROTECT(K = AS_NUMERIC(K));
  PROTECT(lambda = AS_NUMERIC(lambda));
  PROTECT(mu = AS_NUMERIC(mu));
  PROTECT(psi = AS_NUMERIC(psi));
  PROTECT(rho = AS_NUMERIC(rho));
  PROTECT(times = AS_NUMERIC(times));
  PROTECT(ttypes = AS_INTEGER(ttypes));
  PROTECT(add_par = AS_INTEGER(add_par));

  Rboolean parMat = isMatrix(lambda);
  int nrow = 0;
  int parVecLen = 1;
  int quit = 0;
  int surv = INTEGER(add_par)[0];
  double zeroTol = 1e-12;

  /* Perform some checks */
  int vf = (LENGTH(add_par) > 1) ? INTEGER(add_par)[1] : 0;
  if (vf) Rprintf("Verbose mode is on.\n");

  int rs = (LENGTH(add_par) > 2) ? (INTEGER(add_par)[2]>0) : 1;
  if (rs && vf) Rprintf("Rescaling is on.\n");

  int extant = (LENGTH(add_par) > 3) ? INTEGER(add_par)[3] : 0;
  if (vf) Rprintf("%d lineages extant at the root.\n",extant);

  int est_norm = (LENGTH(add_par) > 4) ? INTEGER(add_par)[4] : 0;
  if (est_norm && vf) Rprintf("Forcing estimation of matrix norm.\n");

  // check dimensions
  if (parMat) {
    if (vf) Rprintf("Checking matrix dimensions: ");
    SEXP dim = GET_DIM(lambda);
    nrow = INTEGER(dim)[0];
    parVecLen = INTEGER(dim)[1];
    if (vf) Rprintf("nrow = %d, ncol = %d\n",nrow,parVecLen);
  } else {
    nrow = LENGTH(lambda);
    if (vf) Rprintf("nrow = %d\n",nrow);
  }

  if (LENGTH(K) < parVecLen) quit = 1;

  double maxK = 0.0;
  double* KK = NUMERIC_POINTER(K);

  for (int i = 0; i < parVecLen; ++i) {
    if (maxK < KK[i]) maxK = KK[i];
    if (KK[i] <= 0.0) quit++;
  }

  if (parMat) {
    SEXP dim = GET_DIM(mu);
    if (INTEGER(dim)[0] != nrow) quit++;
    dim = GET_DIM(psi);
    if (INTEGER(dim)[0] != nrow) quit++;
  }

  int maxN = (int) ceil(maxK);
  if (vf) Rprintf("Maximum dimension: N = %d\n",maxN);

  if (nrow != maxN+1) {
    Rprintf("User-supplied functions are not the correct length!");
    Rprintf(" nrow = %d, ",nrow);
    Rprintf(" N = %d, ",maxN+1);
    quit++;
  }

  if (extant < 0 || extant > maxN) {
    Rprintf("Bad root lineages: %d\n",extant);
    quit++;
  }

  if (quit) {
    // not enough columns
    if (vf) Rprintf("Exiting prematurely.\n");
    SEXP p;
    PROTECT(p = NEW_NUMERIC(1));
    REAL(p)[0] = R_NegInf;
    UNPROTECT(9);
    return p;
  }

  double* plambda = NUMERIC_POINTER(lambda);
  double* pmu = NUMERIC_POINTER(mu);
  double* ppsi = NUMERIC_POINTER(psi);
  double* prho = NUMERIC_POINTER(rho);

  double* ptimes  = NUMERIC_POINTER(times);
  int*    pttypes = INTEGER_POINTER(ttypes);
  int     nt      = LENGTH(times);

  int maxExtant = 0;
  int ki = 0;

  SEXP p;
  PROTECT(p = NEW_NUMERIC(maxN+1));  /* +1 for p(0) */
  double* p0 = NUMERIC_POINTER(p);

  double t0 = 0.0;
  double scale = 0.0;

  int i;
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

  int info = -1;

  if (surv) {
    p0[0] = 1.0;
    int i;
    for (i = 1; i <= maxN; ++i) {
      p0[i] = (*prho <= 0.0 || *prho >= 1.0) ? 0.0 : pow(1.-*prho,i);
    }
    expotree_ext(KK,plambda,pmu,ppsi,
                 &parVecLen,&ki,&nt, 
                 ptimes,pttypes,p0,&t0, 
                 &info,&est_norm,&vf,&rs);
    for (i = 0; i <= ceil(KK[parVecLen-1]); ++i) {
      if (p0[i] >= zeroTol) {
        Rprintf("Log survival probability is non-negative! p(%d) = %g\n",i,p0[i]);
      }
      p0[i] = (p0[i] < 0.0) ? log(1.-exp(p0[i])) : R_NegInf;
    }
  } else {
    // set initial value of p
    if (extant == 0) {
      p0[0] = 0.0;
      for (i = 1; i <= maxN; ++i) p0[i] = ppsi[i];
      ki = 1;
      t0 = ptimes[0];
      ptimes = ptimes+1;
      pttypes = pttypes+1;
      --nt;
    } else {
      ki = extant;
      p0[0] = 0.0;
      scale = extant*log(*prho);
      for (i = 1; i <= maxN; ++i) {
        if (i < extant) p0[i] = 0.0;
        else p0[i] = pow(1.-*prho,i-extant);
      }
    }
    expotree_ext(KK,plambda,pmu,ppsi,
                 &parVecLen,&ki,&nt, 
                 ptimes,pttypes,p0,&t0, 
                 &info,&est_norm,&vf,&rs);
    if (info > 0) {
      for (i = 0; i <= ceil(KK[parVecLen-1]); ++i) {
        p0[i] += scale;
      }
    }
  }

  UNPROTECT(10);
  return p;
}

