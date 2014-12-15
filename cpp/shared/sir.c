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

#include "sir.h"

/****************************************************************************/

int sir_index_N(int I, int R, int N) {
  if (I >= 0 && R >= 0 && I+R <= N && I <= N && R <= N) {
    return R*(N+1) - R*(R-1)/2 + I;
  } else {
    return -1;
  }
}

int sir_index(int I, int R, const expo_type* expo) {
  return sir_index_N(I,R,expo->N_max);
}

/****************************************************************************/

int sir_init_index(expo_type* expo) {
  int m = 0;
  if (expo->lookup != NULL) {
    int I, R;
    for (R = 0; R <= expo->N_max; ++R) {
      for (I = 0; I <= expo->N_max-R; ++I) {
        expo->lookup[m] = I;
        expo->lookup[m+expo->dim] = R;
        ++m;
      }
    }
  } else {
    Rprintf("Lookup table not initialized!\n");
  }
  return m;
}

/****************************************************************************/

void sir_start(expo_type* expo, int N) {
  expo->N_max = N;
  expo->dim = sir_index_N(0,N,N)+1;
  expo->num_states = 2;
}

void sir_init_all(expo_type* expo) {
  if (expo->lookup != NULL 
      && expo->lambdaVec != NULL 
      && expo->mat != NULL) 
  {
    sir_init_index(expo);
    sir_init_lambda(expo);
    expo->shift = sir_trace(expo)/expo->dim;
    sir_init_mat(expo);
  } else {
    Rprintf("Allocated memory before initializing!\n");
  }
}

/****************************************************************************/

void sir_init_lambda(expo_type* expo) {
  double S;
  int I, R;
  int m = 0;
  for (R = 0; R <= expo->N_max; ++R) {
    for (I = 0; I <= expo->N_max-R; ++I) {
      S = expo->K - R - I;
      if (S < 0.0) S = 0.0;
      expo->lambdaVec[m] = expo->beta*S/expo->N_max;
      ++m;
    }
  }
}

/****************************************************************************/

void sir_init_mat(expo_type* expo) {
  int I, R;
  int m = 0;

  memset(expo->mat,0,3*expo->dim*sizeof(double));
  memset(expo->mat_i,0,3*expo->dim*sizeof(int));

  if (expo->lookup != NULL) {
    for (m = 0; m < expo->dim; ++m) {
      I = expo->lookup[m];
      R = expo->lookup[m+expo->dim];

      if (I >= expo->ki) {
        expo->mat[m]   = -I*(expo->lambdaVec[m]+expo->mu+expo->psi) - expo->shift;
        expo->mat_i[m] = m;
      }

      if (I > expo->ki && R < expo->N_max) {
        expo->mat[expo->dim+m] = (I-expo->ki)*expo->mu;
        expo->mat_i[expo->dim+m] = sir_index(I-1,R+1,expo);
      }

      if (I >= expo->ki && I+R < expo->N_max) {
        expo->mat[2*expo->dim+m] = (I+expo->ki)*expo->lambdaVec[m];
        expo->mat_i[2*expo->dim+m] = sir_index(I+1,R,expo);
      }
    }
  } else {
    Rprintf("Lookup table not initialized!\n");
  }
  expo->mat_n = m;
}

/****************************************************************************/

void sir_matfunc(char trans, int n1, int n2, 
                 double alpha, double* pin, double* pout,
                 void* pars) 
{
  expo_type* expo = (expo_type*) pars;
  // int dim = sir_index(0,expo->N_max,expo)+1;
  int dim = expo->dim;
  int m;
  int a, b;
  memset(pout,0,expo->dim*sizeof(double));

  int I, R;

  for (m = 0; m < expo->dim; ++m) {
    I = expo->lookup[m];
    R = expo->lookup[m+dim];
    switch (trans) {
      case 'n':
      case 'N':
      default:
        a = sir_index(I-1,R+1,expo);
        b = sir_index(I+1, R ,expo);
        pout[m]  = expo->mat[m] * pin[m];
        if (a >= 0) pout[m] += expo->mat[  dim+m] * pin[a];
        if (b >= 0) pout[m] += expo->mat[2*dim+m] * pin[b];
        pout[m] *= alpha;
        break;
      case 't':
      case 'T':
        a = sir_index(I+1,R-1,expo);
        b = sir_index(I-1, R ,expo);
        pout[m]  = expo->mat[m] * pin[m];
        if (a >= 0) pout[m] += expo->mat[  dim+a] * pin[a];
        if (b >= 0) pout[m] += expo->mat[2*dim+b] * pin[b];
        pout[m] *= alpha;
        break;
    }
  }
}

/****************************************************************************/

double sir_trace(void* pars) {
  expo_type* expo = (expo_type*) pars;
  double trace = 0.0;
  int m = 0;
  int I, R;
  for (R = 0; R <= expo->N_max; ++R) {
    for (I = expo->ki; I <= expo->N_max-R; ++I) {
      // m = sir_index(I,R,expo);
      trace += -I*(expo->lambdaVec[m]+expo->mu+expo->psi);
    }
  }
  return trace;
}

/****************************************************************************/

double sir_inf_norm(void* pars) {
  expo_type* expo = (expo_type*) pars;
  int m = 0;
  double rsum = 0.0;
  double maxrsum = 0.0;
  int dim = sir_index(0,expo->N_max,expo)+1;
  for (m = 0; m < expo->mat_n; ++m) {
    rsum  = fabs(expo->mat[m]);
    rsum += fabs(expo->mat[dim+m]);
    rsum += fabs(expo->mat[2*dim+m]);
    if (rsum > maxrsum) maxrsum = rsum;
  }
  return maxrsum;
}

/****************************************************************************/

double sir_one_norm(void* pars) {
  expo_type* expo = (expo_type*) pars;
  int dim = sir_index(0,expo->N_max,expo)+1;
  double* csum = (double*) calloc(dim,sizeof(double));
  int m = 0;
  for (m = 0; m < expo->mat_n; ++m) {
    csum[expo->mat_i[m]] += fabs(expo->mat[m]);
    csum[expo->mat_i[dim+m]] += fabs(expo->mat[dim+m]);
    csum[expo->mat_i[2*dim+m]] += fabs(expo->mat[2*dim+m]);
  }
  double maxcsum = 0.0;
  for (m = 0; m < dim; ++m) {
    if (csum[m] > maxcsum) maxcsum = csum[m];
  }
  free(csum);
  return maxcsum;
}

/****************************************************************************/

void sir_trans(double* pin, double* pout, expo_type* expo) {
  int I, R;
  int m = 0;
  int n;
  for (R = 0; R <= expo->N_max; ++R) {
    for (I = 0; I < expo->N_max-R; ++I) {
      n = sir_index(I+1,R,expo);
      m = sir_index(I,R,expo);
      pout[m] = 2.0*expo->lambdaVec[m]*pin[n];
      ++m;
    }
    pout[m++] = 0.0; /* I = N-R */
  }
}

/************************************************************/ 

void sir_sample(double* pin, double* pout, expo_type* expo) {
  int I, R;
  int m = 0;
  int n;
  for (R = 0; R < expo->N_max; ++R) {
    pout[m++] = 0.0; /* I = 0 */
    for (I = 1; I <= expo->N_max-R; ++I) {
      n = sir_index(I-1,R+1,expo);
      pout[m] = expo->psi*pin[n];
      ++m;
    }
  }
  /* R = N, I = 0 */
  pout[m++] = 0.0;
}

/************************************************************/ 

void sir_find_nonzero(double* p, expo_type* expo) {
  int I, R;
  int m = 0;
  for (R = 0; R <= expo->N_max; ++R) {
    for (I = 0; I <= expo->N_max-R; ++I) {
      if (p[m] != 0.0) printf("p[%d,%d] = %g\n",I,R,p[m]);
      ++m;
    }
  }
}


