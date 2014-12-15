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

#include "sir_expotree.h"

/************************************************************/ 

void sir_expotree(double* RN, int* Rki, double* Rbeta, double* Rmu,
    double* Rpsi, int* Rn, int* parVecLen, double* times, 
    int* ttypes, double* p, double* t0, int* info, 
    int* estimateNorm, int* Rvflag, int* Rrescale)
{
  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));

  // allocate parameter structure and initialize
  init_expo(expo);

  expo->N       = (int) ceil(RN[0]); /* total population size */
  expo->K       = RN[0];             /* carrying capacity */
  expo->ki      = *Rki;              /* current number of lineages */
  expo->beta    = Rbeta[0];          /* branching rate */
  expo->mu      = Rmu[0];            /* extinction rate */
  expo->psi     = Rpsi[0];           /* sampling rate */
  expo->rescale = *Rrescale;         /* rescale probability vector */
  expo->vflag   = *Rvflag;           /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = *parVecLen;      /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = RN;                /* N parameters */
  expo->betaVec = Rbeta;             /* beta parameters */
  expo->muVec   = Rmu;               /* mu parameters */
  expo->psiVec  = Rpsi;              /* psi parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */
  expo->est_norm = *estimateNorm;    /* force estimation of matrix norm */

  int n = *Rn;                       /* number of events in the tree */

  /* get maximum carrying capacity */
  expo->N_max = max_pop_size(expo);

  /* get dimension of the system */
  expo->dim = sir_index(0,expo->N_max,expo)+1;

  expo->mat = (double*) malloc(3*expo->dim*sizeof(double));
  expo->mat_i = (int*) malloc(3*expo->dim*sizeof(int));

  /* the lookup table is actully only needed when nstates > 1 */
  expo->num_states = 2;
  expo->lookup = (int*) malloc(expo->num_states*expo->dim*sizeof(int));

  /* allocate workspace for lambda */
  expo->lambdaVec = (double*) malloc(expo->dim*sizeof(double));
  expo->psiFun = (double*) malloc(expo->dim*sizeof(double));
  expo->muFun  = (double*) malloc(expo->dim*sizeof(double));

  expo->init_all = &sir_init_all;
  expo->matvec = &sir_matfunc;
  expo->ft = &sir_trans;
  expo->fs = &sir_sample;
  // expo->norm = &sir_one_norm;
  expo->trace = &sir_trace;

  sir_init_index(expo);
  sir_init_all(expo);
  init_mupsi_const(expo);
  init_wrk(expo);

  /* set present time */
  expo->wrk[0] = *t0;

  /* allocate workspace for history */
  expo->p_save = NULL;

  /* call algorithm */
  *info = expoTree(n,times,ttypes,p,expo);

  /* clean up */
  free(expo->mat); expo->mat = NULL;
  free(expo->lambdaVec); expo->lambdaVec = NULL;

  free_wrk(expo);
  free(expo);
  expo = NULL;
}


