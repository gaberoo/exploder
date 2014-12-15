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

#include "expomv.h"

/************************************************************/ 

void rExpoTree(double* RN, int* Rki, double* Rbeta, double* Rmu,
    double* Rpsi, int* Rn, int* parVecLen, double* times, 
    int* ttypes, double* p, double* t0, int* info, 
    int* estimateNorm, int* Rvflag, int* Rrescale) 
{
  /* allocate parameter structure and initialize */
  expo_type* expo = expo_type_alloc();

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

  /* SIS */
  expo_start(expo,max_pop_size(expo));
  expo->init_all = &init_all;
  expo->matvec = &matFuncExpmv;
  
  expo_alloc_all(expo);
  expo->init_all(expo);
  init_mupsi_const(expo);
  init_wrk(expo);

  /* set present time */
  expo->wrk[0] = *t0;

  /* allocate workspace for history */
  // expo->p_save = (double*) malloc(n*expo->dim*sizeof(double));

  /* verbose output */
  if (expo->vflag > 0) {
    Rprintf("Running expoTree with parameters:\n");
    Rprintf(" n    = %d\n",n);
    Rprintf(" ki   = %d\n",expo->ki);
    Rprintf(" Nmax = %d\n",expo->dim);
    int i;
    Rprintf("   K    =");
    for (i = 0; i < expo->parVecLen; ++i) Rprintf(" %8.2f",expo->NVec[i]);
    Rprintf("\n");
    if (! expo->user_funcs) {
      Rprintf("   beta =");
      for (i = 0; i < expo->parVecLen; ++i) Rprintf(" %8.2f",expo->betaVec[i]);
      Rprintf("\n   mu   =");
      for (i = 0; i < expo->parVecLen; ++i) Rprintf(" %8.2f",expo->muVec[i]);
      Rprintf("\n   psi  =");
      for (i = 0; i < expo->parVecLen; ++i) Rprintf(" %8.2f",expo->psiVec[i]);
      Rprintf("\n");
    } else {
      Rprintf("Running with user-supplied functions.\n");
    }
  }

  /* call algorithm */
  *info = expoTree(n,times,ttypes,p,expo);

  /* clean up */
  // free(expo->p_save);

  expo_free_all(expo);
  free_wrk(expo);
  expo_type_free(expo);
}


