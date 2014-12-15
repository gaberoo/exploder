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
  function [c,mv] = normAm(A,m)
  %NORMAM   Estimate of 1-norm of power of matrix.
  %   NORMAM(A,m) estimates norm(A^m,1).
  %   If A has nonnegative elements the estimate is exact.
  %   [C,MV] = NORMAM(A,m) returns the estimate C and the number MV of
  %   matrix-vector products computed involving A or A^*.

  %   Reference: A. H. Al-Mohy and N. J. Higham, A New Scaling and Squaring
  %   Algorithm for the Matrix Exponential, SIAM J. Matrix Anal. Appl. 31(3):
  %   970-989, 2009.

  %   Awad H. Al-Mohy and Nicholas J. Higham, September 7, 2010.
 *
 *
 * Calculate power of matrix-vector product
 *   trans : 'N' no, 'T' transpose
 *   fmv   : matrix-vector product
 *   n     : dimension of square matrix
 *   m     : power
 *   x     : input vector (dim = n)
 *   wrk   : workspace (dim = 2*n)
 *           on output the first n element of wrk are the result of y = (A^m)*x
 */

#include <stdlib.h>
#include <math.h>

#include "../shared/expmv.h"
#include "expocl.h"

void cl_afun_power(char trans, double alpha, int m, clReal* x, 
                   expo_type* expo, expo_cl* ecl)
{
  int i;
  expocl_copy_x(ecl,x);
  for (i = 0; i < m; ++i) {
    ecl->calc(ecl,expo,trans,alpha);
    expocl_copy_xy(ecl);
    // expocl_swap_xy(ecl);
  }
  expocl_finish(ecl);
  expocl_read_x(ecl,x);
}

int expocl_normam(double alpha, int m, double* c, int* mv,
                  double* dwrk, int* iwrk, expo_type* expo, expo_cl* ecl)
{
  int N = expo->dim;

  // dim(wrk) = n + tcol
  // dim(iwrk) = 2*n
  int tcol = 1; // Number of columns used by DLACN1
  char trans = 'N';

  double* mvwrk = dwrk;
  double* v     = mvwrk + 2*N*tcol;
  double* x     = v + N;
  double* xold  = x + N*tcol;
  double* wrk  = xold + N*tcol;
  double* h     = dwrk + tcol;

  int* ind   = iwrk;
  int* indh  = ind + N;

  int kase = 0;
  int info = 0;

  // int iseed[] = { 153, 1673, 2, 3567 };
  int isave[] = { 0, 0, 0 };

  *mv = 0;

  *c = 0.0;

  // dlacn1_(&N,&tcol,v,x,&N,xold,&N,wrk,h,ind,indh,c,&kase,iseed,&info);
  dlacn2_(&N,v,x,iwrk,c,&kase,isave);

  while (kase != 0) {
    if (kase == 1) trans = 'N';
    else if (kase == 2) trans = 'T';
    // call matrix-matrix product
#ifndef DOUBLE
    clReal* b = (clReal*) mvwrk;
    double_to_clreal(N,x,b);
    cl_afun_power(trans,alpha,m,b,expo,ecl);
    clreal_to_double(N,b,x);
#else
    cl_afun_power(trans,alpha,m,x,expo,ecl);
#endif
    // dlacn1_(&N,&tcol,v,x,&N,xold,&N,wrk,h,ind,indh,c,&kase,iseed,&info);
    dlacn2_(&N,v,x,iwrk,c,&kase,isave);
    *mv += m*tcol;
  }

  return 0;
}
