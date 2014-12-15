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

/*
 * Select degree of Taylor approximation. 
 * ======================================
 * 
 * Adapted from MATLAB code of 
 *    Awad H. Al-Mohy and Nicholas J. Higham, October 26, 2010.
 *  Reference: A. H. Al-Mohy and N. J. Higham, Computing the action of
 *  the matrix exponential, with an application to exponential
 *  integrators. MIMS EPrint 2010.30, The University of Manchester, 2010.
 *
 * Gabriel E Leventhal, November 13, 2011.
 *
 * PARAMETERS
 *
 * n (input)    : dimension of matrix
 * fmv (intput) : matrix-matrix product function
 * b (intput)   : matrix b in x = A*b
 * m (intput)   : number of columns of b
 * m_max (input) : default = 8
 * p_max (input) : default = 55
 * prec (input)  : desired precision ('d' = double, 's' = single, 'h' = half)
 * alpha (output) : array of dimension p_max-1
 * eta (output)   : array of dimension p_max
 * M (output)     : approximation matrix (dimension (m_max,p_max))
 * shift (input) : set to true if a shift of the matrix is desired
 *                 (not implemented; shifting on matrix-vector product?)
 * bal (input)   : set to true is balancing of the matrix is desired
 *                 (not implemented; balancing on matrix-vector product?)
 * force_estm (intput) : default = false
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../shared/expmv.h"

#ifdef THETA_HALF
#include "../shared/theta_taylor_half.h"
#elif defined THETA_SINGLE
#include "../shared/theta_taylor_single.h"
#else
#include "../shared/theta_taylor.h"
#endif

#include "expocl.h"

double max_dt(expo_type* expo) {
  /* calculate maximum matrix norm */
  double normLim = 4.0 * theta[expo->m_max-1] * expo->p_max
                   * (expo->p_max+3.0) / (1.0*expo->m_max);
  double normA = expo->norm(expo);
  return normLim/normA;
}

int cl_stdeg(double dt, clReal* b, expo_type* expo, expo_cl* ecl, 
             int force_estm, double* alpha, double* eta)
{
  int m = 1;

  /* check if p and m are acceptable limits */
  if (expo->p_max < 2 
      || expo->m_max + 1 < expo->p_max*(expo->p_max-1) 
      || expo->m_max > 100) 
  {
    memset(expo->tm,0,expo->p_max*expo->m_max*sizeof(double));
    return -1;
  }

  int i;
  int j;

  /* calculate maximum matrix norm */
  double normLim = 4 * theta[expo->m_max-1] * expo->p_max
                   * (expo->p_max+3) / (expo->m_max*m);
  int mv = 0;

  double normA = 0.0;
  double c = 0.0;

  /* do not force reestimation */
  if (!force_estm) {
    normA = dt*expo->norm(expo);
    if (normA > normLim) {
      if (expo->vflag > 1) Rprintf("## ");
      // Rprintf("Matrix norm too large. Decrease time step!\n");
    } else {
      c = normA;
      for (i = 0; i < expo->p_max-1; ++i) alpha[i] = c;
    }
  }

  /* check if norm is acceptable */
  if (force_estm || normA > normLim) {
#ifdef DEBUG
    if (normA > normLim) Rprintf("Matrix norm too large! |A| = %g.",normA);
    // Rprintf("Estimating 1-norm...");
#endif

    if (expo->wrk == NULL || expo->iwrk == NULL || eta == NULL) {
      Rprintf("Must supply workspace to estimate 1-norm.\n");
      return -10;
    }

    /* copy probbility vector from device to the host */
    expocl_read_b(ecl,b);

    int k = 0;
    for (i = 1; i <= expo->p_max; ++i) {
      expocl_normam(dt,i+1,&c,&k,expo->wrk+expo->p_max,expo->iwrk,expo,ecl);
      c = pow(c,1./(i+1.));
      mv += k;
      eta[i] = c;
    }

    for (i = 0; i < expo->p_max-1; ++i) {
      alpha[i] = (eta[i] > eta[i+1]) ? eta[i] : eta[i+1];
    }

#ifdef DEBUG
    // Rprintf("done.\n");
#endif
  }

  for (j = 1; j < expo->p_max; ++j) {
    memset(expo->tm+(j-1)*expo->m_max,0,((j+1)*j-2)*sizeof(double));
    for (i = (j+1)*j-2; i < expo->m_max; ++i) {
      expo->tm[(j-1)*expo->m_max+i] = alpha[j-1]/theta[i];
    }
  }

  return mv;
}

