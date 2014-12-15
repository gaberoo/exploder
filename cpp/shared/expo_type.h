/*
 * Copyright (c) 2012-2013, Gabriel Leventhal, ETH Zurich
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of the <organization> nor the
 *     names of its contributors may be used to endorse or promote products
 *     derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __EXPO_TYPE_H__
#define __EXPO_TYPE_H__

typedef struct _expo {
  int N;        /* maximum population size */
  int ki;       /* current number of lineages */

  int dim;       /* dimension of system */
  int parVecLen; /* number of parameter sets */
  int curPar;    /* current parameter set */

  int vflag;      /* verbosity flag */
  int rescale;    /* rescale probability vector at each iteration */

  int SImodel;    /* use desity-dependent model */

  int offset;
  int useLog;
  int includeZero;

  int p_max;
  int m_max;
  int N_max;

  int max_wrk_steps; /* maximum number of work steps for each interval */

  int est_norm;      /* perform numerical estimation of matrix norm */
  int user_funcs;    /* use user-supplied functions (muFun,psiFun) */

  double beta;       /* infection/speciation rate */
  double mu;         /* recovery/extinction rate */
  double psi;        /* per-lineage sampling rate */
  double rho;        /* proability of discovery at present */
  double K;          /* real-valued carrying capacity */

  double* NVec;      /* carrying capacity parameters */
  double* betaVec;   /* infection rate parameters */
  double* muVec;     /* recovery rate parameters */
  double* psiVec;    /* sampling rate parameters */

  double gamma;
  double shift;

  double tol;     /* numerical tolerance */
  double cutoff;  /* numerical cutoff for zero */

  double mat_inf_norm;
  double mat_one_norm;

  double prec;
  int    full_term;
  int    check_negs;

  int     num_states;
  int*    lookup; /* lookup table for indices */

  double* mat;    /* matrix entries */
  int*    mat_i;  /* row indices */
  int     mat_n;  /* number of non-zero rows */

  double* lambdaVec;
  double* muFun;
  double* psiFun;

  double (*lambda)(int, struct _expo*); /* infection rate function */

  double (*norm)(void*);
  double (*trace)(void*);

  void (*ft)(double*,double*,struct _expo*);
  void (*fs)(double*,double*,struct _expo*);
  void (*fs2)(double*,double*,struct _expo*);

  void   (*matvec)(char,int,int,double,double*,double*,void*);

  void (*init_all)(struct _expo*);

  double* p_save;

  int wlen;
  int iwlen;
  double* wrk;
  int* iwrk;

  double* tm;
} expo_type;

#endif /* __EXPO_TYPE_H__ */
