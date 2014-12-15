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
 * input
 * -----
 *
 * t         : time
 * A         : matrix of derivatives
 * b         : starting vector
 * M         : taylor degree
 * prec      : required accuracy (double, single, half)
 * shift     : shift matrix (default true)
 * bal       : balance matrix (default false)
 *               N = none
 *               P = permute only
 *               S = scale only
 *               B = scale + permute
 * full_term : default false
 * prnt      : default false
 *
 *
 *
 * output
 * ------
 *  mv  : number of matrix-vector products
 *  mvd : (mvd < mv) number of MV products that were used for norm
 *        estimation
 *
 * ORIGINAL FUNCTION DESCRIPTION (MATLAB CODE)
 * ===========================================
 *
 * function [f,s,m,mv,mvd,unA] = ...
 *  %EXPMV   Matrix exponential times vector or matrix.
 *  %   [F,S,M,MV,MVD] = EXPMV(t,A,B,[],PREC) computes EXPM(t*A)*B without
 *  %   explicitly forming EXPM(t*A). PREC is the required accuracy, 'double',
 *  %   'single' or 'half', and defaults to CLASS(A).
 *  %   A total of MV products with A or A^* are used, of which MVD are
 *  %   for norm estimation.
 *  %   The full syntax is
 *  %     [f,s,m,mv,mvd,unA] = expmv(t,A,b,M,prec,shift,bal,full_term,prnt).
 *  %   unA = 1 if the alpha_p were used instead of norm(A).
 *  %   If repeated invocation of EXPMV is required for several values of t
 *  %   or B, it is recommended to provide M as an external parameter as
 *  %   M = SELECT_TAYLOR_DEGREE(A,m_max,p_max,prec,shift,bal,true).
 *  %   This also allows choosing different m_max and p_max.
 *
 *  %   Reference: A. H. Al-Mohy and N. J. Higham, Computing the action of
 *  %   the matrix exponential, with an application to exponential
 *  %   integrators. MIMS EPrint 2010.30, The University of Manchester, 2010.
 *
 *  %   Awad H. Al-Mohy and Nicholas J. Higham, October 26, 2010.
 *
 *
 *
 */
/*
 * PARAMETRS
 *
 * t (input) : scalar in y = exp(tA)*b
 * fmv (intput) : matrix-matrix product function
 * nf (input)   : function to calculate matrix 1-norm
 * tf (input)   : function to calculate trace
 * b (input)    : matrix b in y = exp(tA)*b
 * M (input)    : length of Taylor truncation 
 *                (will be calculated if M <= 0 or M > n)
 * prec (input) : desired precision
 * shift (input) : should the matrix be shifted?
 * bal (intput)  : should the matrix be balanced?
 *
 */

#ifndef __EXPMV_H__
#define __EXPMV_H__

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <R.h>
#include <time.h>

#include "expo_type.h"

typedef void (*matMat)(char trans, int n1, int n2, double alpha, double* x, double* y, void* pars);
typedef double (*normFunc)(void* pars);
typedef double (*traceFunc)(void* pars);

double inf_norm(int n1, int n2, double* A);
double find_absmax_array(int i1, int i2, const double* x);

int expmv(double t, double* b, int recalcm, expo_type* expo);
int select_taylor_degree(double t, double* b, expo_type* expo, int force_estm, 
                         double* alpha, double* eta, double* wrk, int* iwrk);
int normAm(double alpha, int m, double* est, double* dwrk, int* iwrk, expo_type* expo);
void afun_power(char trans, double alpha, int m, double* x, expo_type* expo, double* wrk);
void taylor_length(double tt, const double* tm, const void* expo, double* s, int* tcol);

/* LAPACK ROUTINES */

#if defined(MKL)
#include <mkl.h>
#else
double dlange_(char* norm, int* M, int* N, double* A, int* LDA, void* WORK);
int idamax_(int* n, double* x, int* incx);
void dlacon_(int* n, double* v, double* x, int* isgn, double* est, int* kase);
void dlacn2_(int* n, double* v, double* x, int* isgn, double* est, int* kase, int* isave);
#endif

double dscal_(int* n, double* alpha, double* x, int* incx);
void daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
double dnrm2_(int*,double*,int*);

void dlacn1_(int* n, int* t, double* v, double* x, int* ldx, double* xold, 
    int* ldxold, double* wrk, double* h, int* ind, int* indh, double* est, 
    int* kase, int* iseed, int* info);

#endif
