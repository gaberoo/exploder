/*
 * Copyright (c) 2012-2014, Gabriel Leventhal, ETH Zurich
 * All rights reserved.
 *
 * Original MATLAB implementation:
 * Copyright (c) 2010, Nick Higham and Awad Al-Mohy
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

#ifdef DEBUG
#include <time.h>
clock_t tx = 0;
clock_t t_infnorm = 0;
clock_t t_matvec  = 0;
size_t  c_infnorm = 0;
size_t  c_matvec  = 0;

inline void inc_clock(clock_t* t, clock_t* last_t) {
  *t += clock() - *last_t;
  *last_t = clock();
}
#endif

int cl_expmv(double t, clReal* b, int recalcm,
             double* talpha, expo_type* expo, expo_cl* ecl)
{
  int mv = 0;       /* number of matrix-vector products computed */
  int k = 0;
  double tt;
  int tcol = 1;     /* size of Taylor approximation */

  // get required Taylor truncation (if not set)
  if (recalcm) {
    tt = 1.0;
    mv = cl_stdeg(t, b, expo, ecl, expo->est_norm, talpha, expo->wrk);
    if (mv < 0) {
      Rprintf("Error in stdeg (%d).\n",mv);
      return mv;
    }
  } else {
    tt = t; 
    mv = 0; 
  }

  double s = 1.0;      /* cost per column */

  if (t == 0.0) { 
    tcol = 0; 
  } else {
    taylor_length(tt,expo->tm,expo,&s,&tcol);
  }

  if (tcol == 0) {
    Rprintf("Cannot calculate matrix exponential (under-/overflow?).\n");
    Rprintf("Returned results may be gibberish.\n");
    return -1;
  }

  // Rprintf("= ");

  /* scaling factor eta */
  double eta = (expo->shift != 0.0) ? exp(t*expo->shift/s) : 1.0;

  /* copy b to working vector on device */
  expocl_copy_xb(ecl);

  double c1 = 0.0;
  double c2 = 0.0;
  double bnorm = 0.0;

  int ss;
  for (ss = 1; ss <= s; ++ss) {
    /* get norm of input vector/matrix */
    if (! expo->full_term) c1 = expocl_inf_norm_x(ecl);

    for (k = 1; k <= tcol; ++k) {
      ecl->calc(ecl,expo,'n',t/(s*k)); /* calculate y = A x */
      // expocl_swap_xy(ecl);      
      expocl_copy_xy(ecl);             /* x <=> y */
      ++mv;                            /* increment matrix-vector product counter */
      expocl_add_by(ecl);              /* add x onto b */

     /* check whether convergence has been reached before full truncation */
      if (! expo->full_term) {
        expocl_finish(ecl);             /* make sure the queue has finished */
        c2 = expocl_inf_norm_x(ecl);    /* norm of output vector */
        bnorm = expocl_inf_norm_b(ecl); /* get new infinite norm of the new b */
#ifdef DEBUG
        if (! expo->full_term && expo->vflag > 3) {
          Rprintf("  k=%3d: |b| = %9.2e, c1 = %9.2e, c2 = %9.2e.\n",k,bnorm,c1,c2);
        }
#endif
        /* check if tolerance has been reached already */
        if (c1+c2 <= expo->prec*bnorm) break;
        c1 = c2;
      }
    }

#ifdef DEBUG
    if (! expo->full_term && expo->vflag > 2) {
      Rprintf("m = %3d: %9.2e %9.2e %9.2e %9.2e\n",k,bnorm,c1,c2,(c1+c2)/bnorm);
    }
#endif
    expocl_smult_b(ecl,eta);
    expocl_copy_xb(ecl);
  }

  // expocl_finish(ecl);

  if (expo->iwrk != NULL) {
    expo->iwrk[0] = tcol;
    expo->iwrk[1] = k;
    expo->iwrk[2] = ceil(s);
  }

  return mv;
}

/****************************************************************************/


