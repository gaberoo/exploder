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
#include "expmv.h"

void afun_power(char trans, double alpha, int m, double* x, expo_type* expo, double* wrk)
{
  int i;
  double* y = wrk;
  double* yout = wrk+expo->dim;
  double* ytmp;
  memcpy(y,x,expo->dim*sizeof(double));
  for (i = 0; i < m; ++i) {
    expo->matvec(trans,expo->dim,1,alpha,y,yout,expo);
    // swap pointers
    ytmp = y; 
    y = yout; 
    yout = ytmp;
  }
  memcpy(x,y,expo->dim*sizeof(double));
}


int normAm(double alpha, int m, double* est, double* dwrk, int* iwrk, expo_type* expo)
{
  int tcol = 1; // Number of columns used by DLACN1
  char trans = 'N';

  double* mvwrk = dwrk;
  double* v     = mvwrk + 2*expo->dim;
  double* x     = v + expo->dim;
  double* xold  = x + expo->dim;
  double* wrk   = xold + expo->dim;
  double* h     = wrk + tcol;

  int* ind   = iwrk;
  int* indh  = ind + expo->dim;

  int kase = 0;
  int info = 0;

  int iseed[] = { 153, 1673, 2, 3567 };
  int isave[] = { 0, 0, 0 };

  int mv = 0;

  *est = 0.0;
  // dlacn1_(&expo->dim,&tcol,v,x,&expo->dim,xold,&expo->dim,wrk,h,ind,indh,est,&kase,iseed,&info);
  dlacn2_(&expo->dim,v,x,iwrk,est,&kase,isave);

  while (kase != 0) {
    if (kase == 1) trans = 'N';
    else if (kase == 2) trans = 'T';
    afun_power(trans,alpha,m,x,expo,mvwrk);
    // dlacn1_(&expo->dim,&tcol,v,x,&expo->dim,xold,&expo->dim,wrk,h,ind,indh,est,&kase,iseed,&info);
    dlacn2_(&expo->dim,v,x,iwrk,est,&kase,isave);
    mv += m*tcol;
  }

  return mv;
}
