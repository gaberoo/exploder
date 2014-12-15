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

#include <stdlib.h>
#include "expomv.h"

int integrate_interval(double dt, expo_type* expo, double* b, 
                       double* scale, double* nrm) 
{
  int N = expo->dim;
  int one = 1;

  double* bak = (double*) malloc(expo->dim*sizeof(double));
  memcpy(bak,b,N*sizeof(double));

  int info = 0;

  double wrk_dt = 0.0;
  int wrk_steps = 1;
  int wrk_info = 0;

  double wrk_scale = 0.0;
  double wrk_nrm = 0.0;

  /* TODO: is it really necessary to init_all ? Or only the shift and the
    * matrix ? */
  expo->init_all(expo);

  while (wrk_info == 0) {
    if (wrk_steps > expo->max_wrk_steps) {
      Rprintf("Maximum time step intervals reached.\n");
//#ifdef DEBUG
//      int j;
//      if (! expo->user_funcs) {
//        for (j = 0; j < N; ++j) b[j] = R_NegInf;
//        for (j = 0; j < expo->parVecLen; ++j) Rprintf(" -N %g",expo->NVec[j]);
//        for (j = 0; j < expo->parVecLen; ++j) Rprintf(" -b %g",expo->betaVec[j]);
//        for (j = 0; j < expo->parVecLen; ++j) Rprintf(" -u %g",expo->muVec[j]);
//        for (j = 0; j < expo->parVecLen; ++j) Rprintf(" -s %g",expo->psiVec[j]);
//        for (j = 0; j < N; ++j) if (ttypes[j] == 20) Rprintf(" -S %g",times[j]);
//        Rprintf("\n");
//      }
//#endif
      free(bak);
      return -2;
    }

    wrk_dt = (wrk_steps > 1) ? dt/wrk_steps : dt;
    wrk_scale = 0.0;
    wrk_nrm = 0.0;

    int k;
    for (k = 0; k < wrk_steps; ++k) {
      info = expmv(wrk_dt,b,1,expo);

      // Error during calculation. Return INF.
      if (info < 0) {
        // Rprintf("Error in 'expmv': info = %d\n",info);
        free(bak);
        return -3;
      }

      // check for negative values
      int found_neg = 0;
      int j;
      for (j = 0; j < expo->dim; ++j) {
        if (b[j] < 0.0 || ! isfinite(b[j])) {
          // printf(" --- m = %d, p = %g\n",j,b[j]);
          ++found_neg;
          b[j] = 0.0;
        }
      }

      if (found_neg && expo->vflag > 1) {
        Rprintf("DEBUG: %d negative probabilities found! info = %d\n",found_neg,info);
//        Rprintf("%6d/%02d %8.4f/%8.4f %5d %10.4e % 11.4e | %2d %2d\n",
//                i,ttypes[i],dt,times[i],expo->ki,nrm,scale+wrk_scale,iwrk[0],iwrk[1]);
        *nrm = dnrm2_(&N,b,&one);
        free(bak);
        return -10;
      }

      // calculate norm of vector
      wrk_nrm = dnrm2_(&N,b,&one);

      // validate 2-norm of the vector
      if (wrk_nrm < expo->cutoff || isnan(wrk_nrm)) {
        if (expo->vflag > 1) Rprintf("Vector norm invalid. Aborting.\n");
        free(bak);
        return -4;
      }

      if (expo->rescale) {
        wrk_scale += log(wrk_nrm);
        wrk_nrm = 1./wrk_nrm;
        dscal_(&N,&wrk_nrm,b,&one);
      }

      if (((expo->iwrk[0] > 40 && expo->iwrk[0] < expo->iwrk[1]) 
            && wrk_steps < expo->max_wrk_steps) || found_neg > 0) 
      {
        memcpy(b,bak,N*sizeof(double));
        ++wrk_steps;
        wrk_info = 0;
        if (expo->vflag > 2) {
          Rprintf("Decreasing time step (%d).\n",wrk_steps);
        }
        break;
      } else {
        wrk_info = 1;
      }
    }
  }

  if (expo->rescale) *scale += wrk_scale;

  free(bak);
  return info;
}


