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

#include "expoTreePSave.h"

double expoTreePSave(int parVecLen, double* K, 
                     double* beta, double* mu, double* psi, 
                     double rho, vector<double>& times, vector<int>& ttypes, 
                     int* np, int* ldp, double* psave, 
                     int extant, int est_norm, int vflag, int rescale,
                     int nroot)
{
  double fx = -INFINITY;
  int info = 1;

  int goodParams = check_params(parVecLen,K,beta,mu,psi,rho);

  double maxK = 0.0;
  for (int i = 0; i < parVecLen; ++i) if (maxK < K[i]) maxK = K[i];
  int maxN = (int) ceil(maxK);
  if (vflag) printf("Maximum dimension: N = %d\n",maxN);
  if (nroot < 0 || nroot > maxN) goodParams = 0;

  if (! goodParams) {
    if (vflag > 0) {
      fprintf(stderr,"Ilegal parameters. Returning inf.\n");
      fprintf(stderr,"N = %g, beta = %g, mu = %g, psi = %g, rho = %g\n",
              K[0],beta[0],mu[0],psi[0],rho);
      return -INFINITY;
    }
  } 
   /**************************************************************************/

  /* allocate parameter structure and initialize */
  expo_type* expo = expo_type_alloc();

  expo->N       = (int) ceil(K[0]);  /* total population size */
  expo->K       = K[0];              /* carrying capacity */
  expo->ki      = extant;            /* current number of lineages */
  expo->beta    = beta[0];           /* branching rate */
  expo->mu      = mu[0];             /* extinction rate */
  expo->psi     = psi[0];            /* sampling rate */
  expo->rescale = rescale;           /* rescale probability vector */
  expo->vflag   = vflag;             /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = parVecLen;       /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = K;                 /* N parameters */
  expo->betaVec = beta;              /* beta parameters */
  expo->muVec   = mu;                /* mu parameters */
  expo->psiVec  = psi;               /* psi parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */
  expo->est_norm = est_norm;         /* force estimation of matrix norm */

  int n = times.size();              /* number of events in the tree */

  /* get maximum carrying capacity */
  expo->N_max = max_pop_size(expo);
  expo->dim = expo->N_max+1;

  /* allocate workspace for matrix */
  expo->mat = (double*) malloc(3*expo->dim*sizeof(double));

  /* allocate workspace for lambda */
  expo->lambda = &lambdaSI;
  expo->lambdaVec = (double*) malloc(expo->dim*sizeof(double));
  expo->init_all(expo);

  int ncol = 1;                  /* EXPMV parameters */

  /* calculate required memory and allocate */
  int memlen = (2*expo->dim+(expo->p_max-1)*expo->m_max);
  int wlen   = 2*expo->p_max + (6*ncol+3)*expo->dim + ncol + expo->m_max*(expo->p_max-1);
  int iwlen  = 2*expo->dim + 4;

  double* wrk = (double*) malloc((memlen+wlen)*sizeof(double));
  int* iwrk = (int*) malloc(iwlen*sizeof(int));

  double* ptimes = times.data();
  int*    pttypes = ttypes.data();
  vector<double> p0(maxN+1);

  double t0 = 0.0;
  double scale = 0.0;

  if (*np >= n && *ldp == expo->dim) {
    expo->p_save = psave;
  } else {
    cerr << "Incorrect dimensions for 'psave'. "
         << "Correct dimensions are: "
         << "np  = " << n << ", "
         << "ldp = " << expo->dim << endl;
  }

  // set initial value of p
  if (extant == 0) {
    p0[0] = 0.0;
    for (int m = 1; m <= maxN; ++m) p0[m] = psi[0];
    expo->ki = 1;
    t0 = times[0];
    ptimes = times.data()+1;
    pttypes = ttypes.data()+1;
    n = times.size()-1;
    memcpy(expo->p_save,p0.data(),expo->dim*sizeof(double));
    expo->p_save += expo->dim;
  } 
  else {
    expo->ki = extant;
    p0[0] = 0.0;
    scale = extant*log(rho);
    for (int m = 1; m <= maxN; ++m) {
      if (m < extant) p0[m] = 0.0;
      else p0[m] = pow(1.-rho,m-extant);
    }
    ptimes = times.data();
    pttypes = ttypes.data();
    n = times.size();
  }
  
  /* set present time */
  wrk[0] = t0;

  /* call algorithm */
  expoTree(n,ptimes,pttypes,p0.data(),memlen+wlen,wrk,iwlen,iwrk,expo);
  info = iwrk[0];

  /* clean up */
  expo->p_save = NULL;

  free(expo->mat);       expo->mat = NULL;
  free(expo->lambdaVec); expo->lambdaVec = NULL;
  free(iwrk);            iwrk = NULL;
  free(wrk);             wrk = NULL;

  expo_type_free(expo);

  if (info > 0) {
    fx = p0[nroot+1] + scale;
    if (vflag > 0) fprintf(stderr,"ln(p(1,t)) = %20.12e\n",fx);
  } else {
    if (vflag > 0) fprintf(stderr,"rExpoTree returned %d!\n",info);
    return -INFINITY;
  }

  return fx;
}

