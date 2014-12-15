/*
 * Copyright (c) 2012-2013, Gabriel Leventhal, ETH Zurich
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

#ifndef __EXPO_H__
#define __EXPO_H__

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include "expo_type.h"
#include "R.h"

/***************************************************************************/
/* BLAS/LAPACK                                                             */
/***************************************************************************/

#if defined(MKL)
#include <mkl.h>
#else
/* vector norm */
float snrm2_(int* n, float* x, int* incx);
double dnrm2_(int* n, double* x, int* incx);
/* scale vector */
/* tridiagonal matrix norm */
double dlangt_(char* norm, int* n, double* lower, 
               double* diag, double* upper);
/* tridiagonal matrix-vector multiply */
int dlagtm_(char *trans, int *n, int *nrhs, double *alpha, 
            double *dl, double *d__, double *du, double *x, 
            int *ldx, double *beta, double *b, int *ldb);
#endif

double dscal_(int* n, double* alpha, double* x, int* incx);

/***************************************************************************/

expo_type* init_expo(expo_type* e);
void expo_start(expo_type* e, int N);

expo_type* expo_type_alloc();
void expo_type_free(expo_type* e);

void expo_alloc_all(expo_type* expo);
void expo_free_all(expo_type* expo);

/************************************************************/ 

void init_all(expo_type* expo);
void init_mat(expo_type* expo);
void init_lambda(expo_type* expo);
void init_mupsi_const(expo_type* expo);

void init_wrk(expo_type* expo);
void free_wrk(expo_type* expo);

void expo_make_index(expo_type* expo);
double sis_mat_entry(const expo_type* expo, size_t m, char d);

int max_pop_size(const expo_type* expo);

double lambdaSI(int I, expo_type* expo);
double lambdaInf(int I, expo_type* expo);

double matRowSum(int m, expo_type* expo);
double matColSum(int m, expo_type* expo);
double matColSumShifted(int m, expo_type* expo);

double matFuncOneNorm(void* pars);
double matFuncOneNormShifted(void* pars);
double matFuncInfNorm(void* pars);
double matFuncTrace(void* pars);

double matFuncOneNormAnalytical(void* pars);

void transEvent(double* pin, double* pout, expo_type* expo);
void sampleEvent(double* pin, double* pout, expo_type* expo);
void sampleEventNoShift(double* pin, double* pout, expo_type* expo);

void matFuncExpmv(char trans, int n1, int n2, 
                  double alpha, double* pin, double* pout,
                  void* pars);

#endif // __EXPO_H_
