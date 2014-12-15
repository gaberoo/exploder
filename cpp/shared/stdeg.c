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
 * bal (input)   : set to true is balancing of the matrix is desired
 *                 (not implemented; balancing on matrix-vector product?)
 * force_estm (intput) : default = false
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "expmv.h"

#ifdef THETA_HALF
#include "theta_taylor_half.h"
#elif defined THETA_SINGLE
#include "theta_taylor_single.h"
#else
#include "theta_taylor.h"
#endif

int select_taylor_degree(double t, double* b, expo_type* expo,
                         int force_estm, double* alpha, double* eta,
                         double* wrk, int* iwrk)
{
  if (expo->p_max < 2 || expo->m_max + 1 < expo->p_max*(expo->p_max-1) || expo->m_max > 100) {
    return -1;
  }

  int k = 0;
  double normA = 0.0;
  double c = 0.0;
  int i;
  int j;
  int mv = 0;
  int ncol = 1;

  double normLim = 4*theta[expo->m_max-1]*expo->p_max*(expo->p_max+3)/(expo->m_max*ncol);

  if (!force_estm) normA = t*expo->norm(expo);

  if (!force_estm && normA <= normLim) {
    c = normA;
    for (i = 0; i < expo->p_max-1; ++i) alpha[i] = c;
  } else {
#ifdef DEBUG
    if (normA > normLim) Rprintf("Matrix norm too large! |A| = %g.",normA);
    Rprintf("Estimating 1-norm...");
#endif
    for (i = 1; i <= expo->p_max; ++i) {
      mv += normAm(t,i+1,&c,wrk,iwrk,expo);
      c = pow(c,1./(i+1.));
      mv = mv + k;
      eta[i] = c;
    }

    for (i = 0; i < expo->p_max-1; ++i) {
      alpha[i] = (eta[i] > eta[i+1]) ? eta[i] : eta[i+1];
    }
#ifdef DEBUG
    Rprintf("done.\n");
#endif
  }

  for (j = 1; j < expo->p_max; ++j) {
    memset(expo->tm+(j-1)*expo->m_max,0,((j+1)*j-2)*sizeof(double));
    // for (i = 0; i < (j+1)*j-2; ++i) expo->tm[(j-1)*expo->m_max+i] = 0.0;
    for (i = (j+1)*j-2; i < expo->m_max; ++i) {
      expo->tm[(j-1)*expo->m_max+i] = alpha[j-1]/theta[i];
    }
  }

  return 0;
}


