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

#ifndef __CL_EXPOMV_H__
#define __CL_EXPOMV_H__

#include "../shared/expo.h"
#include "../shared/sir.h"
#include "../shared/expomv.h"
#include "expocl.h"

double max_dt(expo_type* expo);

int cl_stdeg(double dt, clReal* b, expo_type* expo, expo_cl* ecl, 
             int force_estm, double* alpha, double* eta);

int cl_expmv(double t, clReal* b, int recalcm,
             double* talpha, expo_type* expo, expo_cl* ecl);

int clExpoTree(int n, double* times, int* ttypes, double* p, 
               double t, expo_type* expo, expo_cl* ecl);

int cl_integrate_interval(double dt, expo_type* expo, expo_cl* ecl,
                          clReal* b, double* scale, double* nrm);

void cl_update_vector(int ttype, expo_type* expo, expo_cl* ecl);

void clrExpoTree(double* RN, int* Rki, double* Rbeta, double* Rmu,
    double* Rpsi, int* Rn, int* parVecLen, double* times, 
    int* ttypes, double* p, double* t0, int* info, 
    int* estimateNorm, int* Rvflag, int* Rrescale);

#endif
