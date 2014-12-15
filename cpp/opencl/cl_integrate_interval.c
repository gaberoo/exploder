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

#include "cl_expomv.h"

int cl_integrate_interval(double dt, expo_type* expo, expo_cl* ecl,
                          clReal* b, double* scale, double* nrm) 
{
  /* TODO: is it really necessary to init_all ? Or only the shift and the
   * matrix ? */
  expo->init_all(expo);
  /* TODO: may be faster to make matrix on the device, rather than copying? */
  expocl_copy_mat(ecl,expo);

  int wrk_steps = 1;
  double wrk_nrm = 0.0;
  int wrk_mv = 0;

//  int N = expo->dim;
//  int one = 1;

  /* maximum time step allowed */
  double mdt = max_dt(expo);
  if (dt > mdt) {
    wrk_steps = (int) ceil(dt/mdt);
    if (expo->vflag > 1) Rprintf("** ");
    // Rprintf("dt too large. splitting into %d work steps. ",wrk_steps);
    // Rprintf("dt = %f, mdt = %f\n",dt,mdt);
  } else {
    Rprintf(".. ");
  }
  
  double* talpha = (double*) malloc(expo->p_max*sizeof(double));

  if (wrk_steps > expo->max_wrk_steps) { return -102; }

  double wrk_dt = (wrk_steps > 1) ? dt/wrk_steps : dt;

  int k = 0;
  for (k = 0; k < wrk_steps; ++k) {
    int mv = cl_expmv(wrk_dt,b,1,talpha,expo,ecl);

    if (mv <= 0) {
      free(talpha);
      return -103;
    }

    // check for negative values
    if (expo->check_negs) {
      int found_neg = 0;
      found_neg = expocl_neg_vals(ecl);
      if (found_neg) expocl_copy_bx(ecl);

      if (found_neg && expo->vflag > 1) {
        Rprintf("DEBUG: %d negative probabilities found! info = %d\n",found_neg,mv);
#ifdef DEBUG
        wrk_nrm = expocl_two_norm_b(ecl);
        *info = -110;
        free(talpha);
        return;
#endif
      }
    }

    if (expo->rescale) {
      // calculate norm of vector
      wrk_nrm = expocl_two_norm_b(ecl);

      // validate 2-norm of the vector
      if (wrk_nrm < expo->cutoff || isnan(wrk_nrm)) {
        if (expo->vflag > 1) Rprintf("Vector norm invalid. Aborting.\n");
        free(talpha);
        return -104;
      }

      *scale += log(wrk_nrm);
      wrk_nrm = 1./wrk_nrm;
      expocl_smult_b(ecl,wrk_nrm); // dscal_(&N,&wrk_nrm,p0,&one);
    }

    wrk_mv += mv;

    expocl_copy_xb(ecl);
    // expocl_finish(ecl);
  }

  *nrm = wrk_nrm;

  free(talpha);

  return wrk_mv;
}


