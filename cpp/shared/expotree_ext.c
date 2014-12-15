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

#include "expotree_ext.h"

void expotree_ext(double* K, double* lambda, double* mu,
                  double* psi, int* parVecLen, int* ki, int* _n, 
                  double* times, int* ttypes, double* p, double* t0, 
                  int* info, int* estNorm, int* vflag, int* rescale)
{
  /* allocate parameter structure and initialize */
  expo_type* expo = expo_type_alloc();

  expo->N       = (int) ceil(K[0]);  /* total population size */
  expo->K       = K[0];              /* carrying capacity */
  expo->ki      = *ki;               /* current number of lineages */
  expo->rescale = *rescale;          /* rescale probability vector */
  expo->vflag   = *vflag;            /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = *parVecLen;      /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = K;                 /* N parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */
  expo->est_norm = *estNorm;         /* force estimation of matrix norm */

  expo->user_funcs = 1;              /* use user-supplied functions */
  expo->trace = &matFuncTrace_ext;   /* matrix trace */
  expo->fs = &sampleEvent_ext;       /* sampling event */
  expo->init_all = &init_mat_ext;    /* only the matrix needs to be init'ed */

  int n = *_n;                       /* number of events in the tree */

  /* get maximum carrying capacity */
  expo->N_max = max_pop_size(expo);
  int N = expo->N_max+1;

  /* allocate workspace for matrix */
  expo->mat = (double*) malloc(3*N*sizeof(double));

  /* allocate workspace for lambda */
  expo->lambdaVec = lambda;
  expo->muFun     = mu;
  expo->psiFun    = psi;

  expo->init_all(expo);
  init_wrk(expo);

  /* set present time */
  expo->wrk[0] = *t0;

  /* allocate workspace for history */
  expo->p_save = (double*) malloc(n*N*sizeof(double));

  /* call algorithm */
  *info = expoTree(n,times,ttypes,p,expo);

  /* clean up */
  free(expo->p_save); expo->p_save = NULL;
  free(expo->mat); expo->mat = NULL;

  free_wrk(expo);
  expo_type_free(expo);
}

/****************************************************************************/

void init_mat_ext(expo_type* expo) {
  int m;
  int N = expo->N_max+1;
  for (m = 0; m < N; ++m) {
    if (m < expo->ki) {
      expo->mat[m] = 0.0;
      expo->mat[m+N] = 0.0;
      expo->mat[m+2*N] = 0.0;
    } else {
      /* lower diagonal */
      expo->mat[m]     = (m > expo->ki) ? (m-expo->ki)*expo->muFun[m] : 0.0;
      /*    diagonal    */
      expo->mat[N+m]   = -m*(expo->lambdaVec[m]+expo->psiFun[m]+expo->muFun[m]) - expo->shift;
      /* upper diagonal */
      expo->mat[2*N+m] = (m < expo->N) ? (m+expo->ki)*expo->lambdaVec[m] : 0.0;
    }
  }
}

/****************************************************************************/

double matFuncTrace_ext(void* pars) {
  expo_type* expo = (expo_type*) pars;
  double trace = 0.0;
  int m;
  for (m = expo->ki; m < expo->N_max+1; ++m) {
    trace += -m*(expo->lambdaVec[m]+expo->psiFun[m]+expo->muFun[m]);
  }
  return trace;
}

/************************************************************/ 

void sampleEvent_ext(double* pin, double* pout, expo_type* expo) {
  int m;
  pout[0] = 0.0;
  for (m = 1; m <= expo->N_max; ++m) {
    pout[m] = expo->psiFun[m]*pin[m-1];
  }
}


