#include <stdlib.h>
#include <stdio.h>
#include "shared/sir.h"
#include "shared/expmv.h"

int main() {
  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);

  expo->N_max = 10;
  expo->K = 10.0;
  expo->beta = 1.0;
  expo->mu = 0.2;
  expo->psi = 0.2;
  expo->ki = 0;
  expo->vflag = 3;

  int dim = sir_index(0,expo->N_max,expo)+1;
  expo->dim = dim;
  printf("dim = %d\n",dim);

  expo->lambdaVec = (double*) malloc(dim*sizeof(double));
  expo->mat = (double*) malloc(3*dim*sizeof(double));
  expo->mat_i = (int*) malloc(3*dim*sizeof(int));

  double* pin  = (double*) calloc(dim,sizeof(double));
  double* pout = (double*) calloc(dim,sizeof(double));

  int m = sir_index(0,10,expo);
  pin[m] = 1.0;
  // for (m = 0; m < sir_index(0,expo->N_max,expo); ++m) pin[m] = 1.0/m;

  expo->init_all = &sir_init_all;
  expo->matvec = &sir_matfunc;
  expo->ft = &sir_trans;
  expo->fs = &sir_sample;
  expo->norm = &sir_one_norm;
  expo->trace = &sir_trace;

  expo->est_norm = 0;
  expo->offset = 0;

  expo->init_all(expo);

  /* calculate required memory and check */
  int ncol   = 1;

  /* calculate required memory and allocate */
  int memlen = (2*expo->dim+(expo->p_max-1)*expo->m_max);
  int wlen   = 2*expo->p_max + (6*ncol+3)*expo->dim + ncol + expo->m_max*(expo->p_max-1);
  int iwlen  = 2*expo->dim + 4;

  double* wrk = (double*) malloc((memlen+wlen)*sizeof(double));
  int* iwrk = (int*) malloc(iwlen*sizeof(int));

  /* set pointers for easier access to variables */
  double* p0 = wrk;
  double* pT = wrk + expo->dim;
  double* tm = pT + expo->dim;
  double* expowrk = wrk + memlen;

  int info = 100;

  double wrk_dt = 0.1;
  memcpy(p0,pin,expo->dim*sizeof(double));

  expmv(wrk_dt, expo->dim, expo->matvec, expo->norm, expo->trace, p0,
        1, expo->m_max, expo->p_max, tm, 1, 'd', expo->shift, 0, 0,
        expo->vflag, &info, expo->est_norm, wlen, expowrk, iwlen,
        iwrk, expo);

  printf("Number of matrix-vector products = %d\n",iwrk[2]);
  printf("info = %d\n",info);
  printf("d[%d] = %g\n",expo->dim-1,expo->mat[expo->dim-1]);

  m = 0;
  int I, R;
  for (R = 0; R <= expo->N_max; ++R) {
    for (I = 0; I <= expo->N_max-R; ++I) {
      printf("p(%2d|%2d,%2d) = % 12g --> % 12g\n",m,I,R,pin[m],p0[m]);
      ++m;
    }
  }

  free(wrk);
  free(iwrk);

  free(pin);
  free(pout);

  free(expo->lambdaVec);
  free(expo->mat_i);
  free(expo->mat);

  return 0;
}


