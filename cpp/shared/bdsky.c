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

#include <stdio.h>

#include "bdsky.h"

void count_lineages(int n, const int* ttypes, int nroot, int* lineages) {
  int i;
  int lins = nroot;
  for (i = 0; i < n; ++i) {
    switch (ttypes[i]) {
      case 0:
      case 2:
      case 4:
      case 10:
        lineages[i] = lins++;
        break;
      case 1:
      case 11:
        lineages[i] = lins--;
        break;
      default:
        lineages[i] = lins;
        break;
    }
  }
}

void find_intervals(int n, const double* x, int m, const double* y, int* ival);

double bdsky_eval(const double* beta, const double* mu, const double* psi, const double* rho,
                  int nshifts, const double* shifts, int n, const double* times, const int* ttypes, 
                  int survival, int vflag) 
{
  // calculate lineages through time at time shifts
  int* lineages = (int*) calloc(n,sizeof(int));
  count_lineages(n,ttypes,0,lineages);

  int* ival = (int*) calloc(n,sizeof(int));
  find_intervals(nshifts,shifts,n,times,ival);

  int* shift_ival = (int*) calloc(nshifts+1,sizeof(int));
  find_intervals(n,times,nshifts,shifts,shift_ival);
  shift_ival[nshifts] = n-1;

  double* A = (double*) malloc((nshifts+1)*sizeof(double));
  double* B = (double*) malloc((nshifts+1)*sizeof(double));
  double* p = (double*) malloc((nshifts+2)*sizeof(double));

  double t0;
  double ti = 0.0;
  p[0] = 1.0;

  int i;

  // recursively calculate Ai, Bi and pi(t[i-1])
  for (i = 0; i <= nshifts; ++i) {
    t0 = (i < nshifts) ? -shifts[i] : -times[n-1];
    A[i] = bdsky_A(beta[i],mu[i],psi[i]);
    B[i] = bdsky_B(beta[i],mu[i],psi[i],rho[i],A[i],p[i]);
    p[i+1] = bdsky_p(beta[i],mu[i],psi[i],rho[i],t0,ti,A[i],B[i]);
    if (vflag) printf("A(%d) = % f, B(%d) = % f, p(%d) = % f\n",i,A[i],i,B[i],i+1,p[i+1]);
    ti = t0;
  }

  double loglik = 0.0;
  double lnq = 0.0;

  int k;
  for (i = 0; i < n; ++i) {
    k = ival[i];
    ti = (k > 0) ? -shifts[k-1] : 0.0;
    lnq = 0.0;

    switch (ttypes[i]) {
      case 1:
        lnq = bdsky_lnq(-times[i],ti,A[k],B[k]);
        loglik += M_LN2 + log(beta[k]) + lnq;
        break;

      case 0:
        lnq = bdsky_lnq(-times[i],ti,A[k],B[k]);
        loglik += log(psi[k]) - lnq;
        break;

      default:
        break;
    }

    if (vflag) printf("ti(%d) = % f, lnq = % f\n",k,ti,lnq);
  }

  if (vflag) printf("loglik = %f\n",loglik);

  ti = 0.0;
  for (k = 0; k <= nshifts; ++k) {
    t0 = (k < nshifts) ? -shifts[k] : -times[n-1];
    lnq = bdsky_lnq(t0,ti,A[k],B[k]);
    if (vflag) {
      printf("lnq(%.1f|%d) = % f, ival = %d, lins = %d\n",t0,k,lnq,shift_ival[k],lineages[shift_ival[k]]);
    }
    loglik += lineages[shift_ival[k]]*lnq;
    ti = t0;
  }

  k = nshifts;
  // loglik -= M_LN2 + log(beta[nshifts]);

  if (vflag) printf("loglik = %f\n",loglik);

  if (survival) {
    k = nshifts;
    t0 = -times[n-1];
    ti = (k > 0) ? -shifts[k-1] : 0.0;
    double p0 = bdsky_p(beta[k],mu[k],psi[k],rho[k],t0,ti,A[k],B[k]);
    if (vflag) printf("t0 = %f, ti = %f, p0 = %f, ll = %f\n",t0,ti,p0,log(1-p0));
    loglik -= log(1.0-p0);
  }

  free(A);
  free(B);
  free(p);

  free(lineages);
  free(shift_ival);
  free(ival);

  return loglik;
}

double bdsky_A(double beta, double mu, double psi) {
  double sum = beta-mu-psi;
  return sqrt(sum*sum + 4*beta*psi);
}

double bdsky_B(double beta, double mu, double psi, double rho, double Ai, double pi_neg1) {
  return ((1-2*(1-rho)*pi_neg1)*beta+mu+psi)/Ai;
}

double bdsky_p(double beta, double mu, double psi, double rho,
               double t0, double ti, double Ai, double Bi) 
{
  double a = beta+mu+psi;
  double eA = exp(Ai*(ti-t0));
  double b = eA*(1.+Bi);
  double num = b-(1.-Bi);
  double den = b+(1.-Bi);
  return (a-Ai*num/den)/(2.*beta);
}

double bdsky_q(double t0, double ti, double Ai, double Bi) {
  double eA = exp(-Ai*(t0-ti));
  double a = (eA*(1.+Bi)+(1.-Bi));
  return 4.0*eA/(a*a);
}

double bdsky_lnq(double t0, double ti, double Ai, double Bi) {
  double ln_eA = Ai*(ti-t0);
  return 2.*M_LN2 + ln_eA - 2.*(ln_eA + log(1.+Bi) + log( 1.+(1.-Bi)/(exp(ln_eA)*(1.+Bi)) ));
}

