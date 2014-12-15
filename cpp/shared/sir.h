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

#ifndef __SIR_H__
#define __SIR_H__

#include "expo.h"

/* Storage format for 2D vector
 * ----------------------------
 *
 * (I,R):
 *
 *   (0,0)   (0,1)   (0,2)   (0,3) ... (0,N-1) (0,N)
 *   (1,0)   (1,1)   (1,2)   (1,3) ... (1,N-1)
 *    ...     ...     ...     ...
 *  (N-2,0) (N-2,1) (N-2,2)
 *  (N-1,0) (N-1,1)
 *   (N,0)
 *
 * Flatten:
 *
 *  (0,0) (1,0) ... (N,0) (0,1) (1,1) ... (N-1,1) (0,2) ... (N-2,2) ...
 *  \_________  ________/ \__________  _________/ \______  _______/
 *            \/                     \/                  \/
 *       N+1 elements            N elements          N-1 elements
 *  
 *  Introduce m = N+1, then (I,R),
 *
 *  R = 0 =>   m elements
 *  R = 1 => m-1 elements
 *  R = 2 => m-2 elements
 *  ...
 *
 *  For example, 
 *
 *  (I,2),
 *    m + (m-1) + I = 2m-1 + I
 *
 *  (I,3),
 *    m + (m-1) + (m-2) + I = 3m-3 + I
 *
 *  (I,4),
 *    m + (m-1) + (m-2) + (m-3) + I = 4m-6 + I
 *
 *  (I,5),
 *    m + (m-1) + (m-2) + (m-3) + (m-4) + I = 5m-10 + I
 * 
 *  (I,R) => Rm - (0+1+2+...+(R-1)) + I = Rm - (R-1)R/2 + I
 *
 *  Thus, 
 *  
 *    (I,R) => R (N+1) - R (R-1)/2 + I
 *
 *  Tests:
 *
 *    ( 0 , 0 ) => 0
 *    ( N , 0 ) => N
 *    ( 0 , 1 ) => N+1 - 0 = N+1
 *    (N-1, 1 ) => N+1 + N-1 = 2N
 *    ( 0 , 2 ) => 2(N+1) - 2(1)/2 + 0 = 2N+2 - 1 = 2N+1
 *    (N-2, 2 ) => 2(N+1) - 2(1)/2 + N-2 = 2N+2-1 + N-2
 *    ( 0 , N ) => N(N+1) - N(N-1)/2 = N^2 + N - N^2/2 + N/2
 *                   = N^2/2 + 3N/2 = N(N+3)/2
 *
 * */

int sir_index(int I, int R, const expo_type* pars);
int sir_index_N(int I, int R, int N);

void sir_start(expo_type* expo, int N);
int sir_init_index(expo_type* expo);

/*
 * p(I,R;t+dt) = (1-(l[I,R]-m-psi)dt) p(I,R;t)
 *               + p(I-1,R+1;t) (I-ki) (m dt)
 *               + p(I+1, R ;t) (I-ki) (l[I,R] dt)
 *               + p(I+1, R ;t)  2 ki  (l[I,R] dt)
 *
 * dp(I,R)/dt =  -(l[I,R]-m-psi)dt) p(I,R;t)
 *               + p(I-1,R+1;t) (I-ki) (m dt)
 *               + p(I+1, R ;t) (I-ki) (l[I,R] dt)
 *               + p(I+1, R ;t)  2 ki  (l[I,R] dt)
 *
 * dp/dt = A p => A(i,j) = ???
 *
 * */

void sir_matfunc(char trans, int n1, int n2, 
                 double alpha, double* pin, double* pout,
                 void* pars);

void sir_init_all(expo_type* expo);
void sir_init_mat(expo_type* expo);
void sir_init_lambda(expo_type* expo);

double sir_trace(void* pars);
double sir_one_norm(void* pars);
double sir_inf_norm(void* pars);

void sir_trans(double* pin, double* pout, expo_type* expo);
void sir_sample(double* pin, double* pout, expo_type* expo);

void sir_find_nonzero(double* p, expo_type* expo);

#endif
