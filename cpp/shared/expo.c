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

#include "expo.h"

expo_type* init_expo(expo_type* e) {
  e->N = 0;
  e->ki = 0; 
  e->beta = 0.0;
  e->mu = 0.0;
  e->psi = 0.0;
  e->rho = 0.0;
  e->K = 0.0;
  e->dim = 0;
  e->parVecLen = 0;
  e->curPar = 0;
  e->NVec = NULL;
  e->betaVec = NULL;
  e->muVec = NULL;
  e->psiVec = NULL;
  e->gamma = 0.0;
  e->shift = 0.0;
  e->vflag = 0;
  e->rescale = 1;
  e->tol = 0.0;
  e->cutoff = 0.0;
  e->SImodel = 1;
  e->lambda = &lambdaSI;
  e->offset = 0;
  e->useLog = 0;
  e->includeZero = 0;
  e->mat_inf_norm = 0.0;
  e->mat_one_norm = 0.0;
  e->mat = NULL;
  e->mat_i = NULL;
  e->mat_n = 0;
  e->num_states = 1;
  e->lookup = NULL;
  e->matvec = &matFuncExpmv;
  e->lambdaVec = NULL;
  e->muFun = NULL;
  e->psiFun = NULL;
  e->norm = &matFuncOneNormAnalytical;
  e->trace = &matFuncTrace;
  e->ft = &transEvent;
  e->fs = &sampleEvent;
  e->fs2 = &sampleEventNoShift;
  e->p_max = 8;  // +1 than used in the loop calculation
  e->m_max = 55; // +1 than used in the loop calculation
  e->N_max = 0;
  e->est_norm = 0;
  e->init_all = &init_all;
  e->p_save = NULL;
  e->init_all = &init_all;
  e->user_funcs = 0;
  e->wlen = 0;
  e->iwlen = 0;
  e->wrk = NULL;
  e->iwrk = NULL;
  e->tm = NULL;
  e->max_wrk_steps = 100;
  e->full_term = 0;
#if defined(HALF)
  e->prec = 9.765625e-4;
#elif defined(DOUBLE)
  e->prec = 1.1102230246251565404236316680908203125e-16; 
#else
  e->prec = 5.9604644775390625e-7;
#endif
  return e;
}

/****************************************************************************/

expo_type* expo_type_alloc() {
  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));
  return init_expo(expo);
}

void expo_type_free(expo_type* expo) {
  free(expo);
  expo = NULL;
}

void expo_make_index(expo_type* expo) {
  int m;
  for (m = 0; m <= expo->N_max; ++m) {
    expo->lookup[m] = m;
  }
}

/************************************************************/ 

void expo_start(expo_type* expo, int N) {
  expo->N_max = N;
  expo->num_states = 1;
  expo->dim = expo->N_max + 1;
}

void expo_alloc_all(expo_type* expo) {
  if (expo->dim < 1) Rprintf("Please initialize 'expo' first.\n");
  else {
    expo->mat = (double*) malloc(3*expo->dim*sizeof(double));
    expo->mat_i = (int*) malloc(3*expo->dim*sizeof(int));

    /* the lookup table is actully only needed when nstates > 1 */
    expo->lookup = (int*) malloc(expo->num_states*expo->dim*sizeof(int));

    /* allocate workspace for lambda */
    expo->lambdaVec = (double*) calloc(expo->dim,sizeof(double));
    expo->psiFun    = (double*) calloc(expo->dim,sizeof(double));
    expo->muFun     = (double*) calloc(expo->dim,sizeof(double));
  }
}

void expo_free_all(expo_type* expo) {
  free(expo->mat);       expo->mat = NULL;
  free(expo->mat_i);     expo->mat_i = NULL;
  free(expo->lookup);    expo->lookup = NULL;
  free(expo->lambdaVec); expo->lambdaVec = NULL;
  free(expo->psiFun);    expo->psiFun = NULL;
  free(expo->muFun);     expo->muFun = NULL;
}

/************************************************************/ 
/* Initialte matrix                                         */
/************************************************************/ 

void init_mat(expo_type* expo) {
  int m;
  int N = expo->N_max+1;
  for (m = 0; m < N; ++m) {
    if (m < expo->ki) {
      expo->mat[m] = 0.0;
      expo->mat[m+N] = 0.0;
      expo->mat[m+2*N] = 0.0;
    } else {
      expo->mat[m]     = sis_mat_entry(expo,m,'l');
      expo->mat[N+m]   = sis_mat_entry(expo,m,'d');
      expo->mat[2*N+m] = sis_mat_entry(expo,m,'u');
    }
  }
}

double sis_mat_entry(const expo_type* expo, size_t m, char d) {
  if (m < expo->ki || m > expo->N_max) return 0.0;

  switch (d) {
    case 'd':
      return -(1.0*m)*(expo->lambdaVec[m]+expo->psi+expo->mu) - expo->shift;
      break;

    case 'l':
      return (m > expo->ki) ? (1.0*m-expo->ki)*expo->mu : 0.0;
      break;
      
    case 'u':
      return (m < expo->N) ? (1.0*m+expo->ki)*expo->lambdaVec[m] : 0.0;
      break;

    default:
      return 0.0;
      break;
  }
}

/************************************************************/ 
/* Initialte lambda vector                                  */
/************************************************************/ 

void init_lambda(expo_type* expo) {
  int m;
  for (m = 0; m <= expo->N_max; ++m) {
    expo->lambdaVec[m] = expo->lambda(m,expo);
  }
}

void init_mupsi_const(expo_type* expo) {
  int m;
  for (m = 0; m < expo->dim; ++m) {
    expo->muFun[m] = expo->mu;
    expo->psiFun[m] = expo->psi;
  }
}

void init_all(expo_type* expo) {
  if (expo->lookup != NULL 
      && expo->lambdaVec != NULL 
      && expo->mat != NULL) 
  {
    expo_make_index(expo);
    init_lambda(expo);
    expo->shift = expo->trace(expo)/expo->dim;  /* shift matrix */
    init_mat(expo);
  } else {
    Rprintf("Allocated memory before initializing!\n");
  }
};

/************************************************************/ 
/* Workspace                                                */
/************************************************************/ 

void init_wrk(expo_type* expo) {
  expo->iwlen  = 2*expo->dim + 4;
  expo->wlen = 2*expo->p_max + 9*expo->dim + 1 + expo->m_max*expo->p_max
               + (2*expo->dim+(expo->p_max-1)*expo->m_max);
  expo->wrk = (double*) malloc(expo->wlen*sizeof(double));
  expo->iwrk = (int*) malloc(expo->iwlen*sizeof(int));
  expo->tm = (double*) calloc(expo->p_max*expo->m_max,sizeof(double));
}

void free_wrk(expo_type* expo) {
  if (expo->wrk != NULL)  { 
    free(expo->wrk); 
    expo->wrk = NULL; 
    expo->wlen = 0;
  }
  if (expo->iwrk != NULL) { 
    free(expo->iwrk); 
    expo->iwrk = NULL; 
    expo->iwlen = 0;
  }
  if (expo->tm != NULL) { 
    free(expo->tm); 
    expo->tm = NULL; 
  }
}

/************************************************************/ 
/* Lambda functions                                         */
/************************************************************/ 

double lambdaSI(int I, expo_type* expo) { 
  if (I > 0) {
    double l = expo->beta*(1.-I/expo->K);
    return (l > 0.0) ? l : 0.0;
  } else {
    return 0.0;
  }
}

/************************************************************/ 

double lambdaInf(int I, expo_type* expo) { 
  if (I > 0 && I <= expo->N) 
    return expo->beta;
  else 
    return 0.0;
}

/************************************************************/ 
/* Column and row sums                                      */
/************************************************************/ 

double matRowSum(int m, expo_type* expo) {
  int N = expo->N_max+1;
  const double *ld = expo->mat;
  const double *d  = ld + N;
  const double *ud = d + N;
  if (m >= 0 && m < N) return ld[m] + d[m] + ud[m];
  else return 0.0;
}

/************************************************************/ 

double matColSum(int m, expo_type* expo) {
  int N = expo->N_max+1;
  const double *ld = expo->mat;
  const double *d  = ld + N;
  const double *ud = d + N;
  double cs = 0.0;
  if (m >= expo->ki && m < N) {
    cs = d[m];
    if (m > 0) cs += ud[m-1];
    if (m < expo->N_max) cs += ld[m+1];
  }
  return cs;
}

/************************************************************/ 
/* Matrix-vector product required for EXPMV                 */
/* - matrix-matrix product                                  */
/* - shifted matrix-matrix product                          */
/* - shifted matrix norm                                    */
/* - matrix trace                                           */
/************************************************************/ 

void matFuncExpmv(char trans, int n1, int n2, 
                  double alpha, double* pin, double* pout,
                  void* pars) 
{
  expo_type* expo = (expo_type*) pars;
  int N = expo->N_max+1;

#if defined(DLAGTM)
  double done  = 1.0;
  double dzero = 0.0;
  int ione = 1;
  dlagtm_(&trans,&N,&ione,&done,
          expo->mat+1,expo->mat+N,expo->mat+2*N,pin, 
          &N,&dzero,pout,&N);
  dscal_(&N,&alpha,pout,&ione);

#elif defined(DDIAGEMV)
  int idiag[] = { -1, 0, 1 };
  int ndiag = 3;
  int ione = 1;
  mkl_ddiagemv(&trans,&N,expo->mat,&N,idiag,&ndiag,pin,pout);
  dscal_(&N,&alpha,pout,&ione);

#else
  double aa = 0.0;
  double bb = 0.0;
  double cc = 0.0;
  int m;
  memset(pout,0,expo->ki*sizeof(double));
  for (m = expo->ki; m < N-1; ++m) {
    switch (trans) {
      case 'n':
      case 'N':
      default:
        cc = expo->mat[m]     * pin[m-1];
        aa = expo->mat[N+m]   * pin[m];
        bb = expo->mat[2*N+m] * pin[m+1];
        break;
      case 't':
      case 'T':
        if (m > 0) bb = expo->mat[2*N+(m-1)]*pin[m-1];
        aa = expo->mat[N+m]*pin[m];
        cc = expo->mat[m+1]*pin[m+1];
        break;
    }
    pout[m] = alpha*(aa+bb+cc);
  }
  switch (trans) {
    case 'n':
    case 'N':
    default:
      cc = expo->mat[m]*pin[m-1];
      aa = expo->mat[N+m]*pin[m];
      break;
    case 't':
    case 'T':
      cc = expo->mat[2*N+(m-1)]*pin[m-1];
      aa = expo->mat[N+m]*pin[m];
      break;
  }
  pout[N-1] = alpha*(aa+cc);
#endif
}

/************************************************************/ 

double matFuncTrace(void* pars) {
  expo_type* expo = (expo_type*) pars;
  double trace = 0.0;
  int m = 0;
  for (m = expo->ki; m < expo->N_max+1; ++m) {
    trace += -m*(expo->lambdaVec[m]+expo->psi+expo->mu);
  }
  return trace;
}

/************************************************************/ 

double matFuncOneNorm(void* pars) {
  expo_type* expo = (expo_type*) pars;
  char norm = 'o';
  int N = expo->N_max+1;
  return dlangt_(&norm,&N,expo->mat+1,expo->mat+N,expo->mat+2*N);
}

/************************************************************/ 

double matFuncInfNorm(void* pars) {
  expo_type* expo = (expo_type*) pars;
  char norm = 'i';
  int N = expo->N_max+1;
  return dlangt_(&norm,&N,expo->mat+1,expo->mat+N,expo->mat+2*N);
}

/************************************************************/ 

double matFuncOneNormAnalytical(void* pars) {
  expo_type* expo = (expo_type*) pars;

  double N = expo->N;
  double m = 0.5*((N-expo->ki+1.)+N/expo->beta*(.5*expo->psi+expo->mu));
  ++m;

  /* slope = zero */
  int m1 = floor(m);
  double c1 = fabs(sis_mat_entry(expo,m1,'d'));
  c1 += fabs(sis_mat_entry(expo,m1-1,'u'));
  c1 += fabs(sis_mat_entry(expo,m1+1,'l'));

  int m2 = floor(m+1);
  double c2 = fabs(sis_mat_entry(expo,m2,'d'));
  c2 += fabs(sis_mat_entry(expo,m2-1,'u'));
  c2 += fabs(sis_mat_entry(expo,m2+1,'l'));

  /* max */
  int m3 = expo->dim-1;
  double c3 = fabs(sis_mat_entry(expo,m3,'d'));
  c3 += fabs(sis_mat_entry(expo,m3-1,'u'));
  c3 += fabs(sis_mat_entry(expo,m3+1,'l'));

  double c4 = fabs(sis_mat_entry(expo,m3-1,'d'));
  c4 += fabs(sis_mat_entry(expo,m3-2,'u'));
  c4 += fabs(sis_mat_entry(expo,m3,'l'));

  /* min */
  int m5 = expo->ki;
  double c5 = fabs(sis_mat_entry(expo,m5,'d'));
  c5 += fabs(sis_mat_entry(expo,m5-1,'u'));
  c5 += fabs(sis_mat_entry(expo,m5+1,'l'));

  double c6 = fabs(sis_mat_entry(expo,m5+1,'d'));
  c6 += fabs(sis_mat_entry(expo,m5,'u'));
  c6 += fabs(sis_mat_entry(expo,m5+2,'l'));

//  printf("m* = %g\n",m);
//  printf("%g %g %g %g\n",c1,c2,c3,c4);

  return fmax(fmax(fmax(c1,c2),fmax(c3,c4)),fmax(c5,c6));
}

/************************************************************/ 

double matFuncInfNormAnalytical(void* pars) {
  expo_type* expo = (expo_type*) pars;
  double maxRow = .5*expo->N + .25*expo->ki 
                  + .25*expo->N/expo->beta*(expo->mu+expo->psi);
  int maxRi = rint(maxRow);
  if (maxRi < expo->ki) maxRi = expo->ki;
  if (maxRi > expo->N) maxRi = expo->N;
  return matRowSum(maxRi,expo);
}

/************************************************************/ 

void transEvent(double* pin, double* pout, expo_type* expo) {
  int m;
  for (m = 0; m < expo->N_max; ++m) {
    pout[m] = 2.0*expo->lambdaVec[m]*pin[m+1];
  }
  pout[expo->N_max] = 0.0;
}

/************************************************************/ 

void sampleEvent(double* pin, double* pout, expo_type* expo) {
  int m;
  pout[0] = 0.0;
  for (m = 1; m <= expo->N_max; ++m) {
    pout[m] = expo->psi*pin[m-1];
  }
  // pout[expo->N_max] = 0.0;
}

/************************************************************/ 

void sampleEventNoShift(double* pin, double* pout, 
    expo_type* expo) {
  int m;

  for (m = 0; m < expo->ki; ++m) 
    pout[m] = 0.0;

  for (m = expo->ki; m <= expo->N_max; ++m) 
    pout[m] = expo->psi*pin[m];
}

/************************************************************/ 

int max_pop_size(const expo_type* expo) {
  /* get maximum carrying capacity */
  int N = (int) ceil(expo->NVec[0]);
  int N2 = 0;
  int j;
  for (j = 0; j < expo->parVecLen; ++j) {
    N2 = (int) ceil(expo->NVec[j]);
    if (N < N2) N = N2;
  }
  return N;
}

/************************************************************/ 

