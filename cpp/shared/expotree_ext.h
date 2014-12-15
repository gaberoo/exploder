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

#ifndef __EXPOTREE_EXT_H__
#define __EXPOTREE_EXT_H__

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>

#include "expo_type.h"
#include "expo.h"
#include "expmv.h"
#include "expomv.h"

/* Run expoTree with custom lambda, mu and psi functions 
 *
 * Parameter vectors:
 *
 *   N        (double[parVecLen])
 *
 * maxN is the smallest integer > 0, s.t. all N <= maxN, i.e. an integer that
 * is larger than all entries in N. This defines the maximum possible
 * population size.
 *
 *   beta     (double[maxN x parVecLen])
 *   mu       (double[maxN x parVecLen])
 *   psi      (double[maxN x parVecLen])
 *
 * The entries of the vectors are such that f[j*maxN+i] is the
 * value of the paramter for _i_ infected individuals in the parameter
 * interval j.
 *
 * parVecLen  (int)          : number of parameter intervals
 * ki         (int)          : number of extant lineages at the start
 * n          (int)          : length of the times vector (including shifts)
 * times      (double[n])    : event times
 * ttypes     (double[n])    : event types
 * p          (double[maxN]) : probability vector
 * t0         (double)       : time at the start
 * info       (int)          : execution information
 * estNorm    (int)          : empirically estimate the matrix norm
 * vflag      (int)          : verbosity level
 * rescale    (int)          : rescale probability vector during run
 *
 * */

void expotree_ext(double* K, double* beta, double* mu,
                  double* psi, int* parVecLen, int* ki, int* n, 
                  double* times, int* ttypes, double* p, double* t0, 
                  int* info, int* estNorm, int* vflag, int* rescale);

void init_mat_ext(expo_type* expo);

void sampleEvent_ext(double* pin, double* pout, expo_type* expo);
double matFuncTrace_ext(void* pars);

#endif
