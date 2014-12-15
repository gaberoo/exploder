#include <stdlib.h>

#include "shared/sir.h"

void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, 
            int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);

int main() {
  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);
  sir_start(expo,5);
  expo->vflag = 3;

  expo->K = 1.0*expo->N_max;
  expo->beta = 1.0;
  expo->mu = 0.2;
  expo->psi = 0.2;
  expo->ki = 1;

  printf("dim = %d\n",expo->dim);
  printf(" N  = %d\n",expo->N_max);

  expo_alloc_all(expo);

  int one = 1;
  int I, R;

  double* pin  = (double*) calloc(expo->dim,sizeof(double));
  double* pout = (double*) calloc(expo->dim,sizeof(double));
  double* pout2 = (double*) calloc(expo->dim,sizeof(double));

  int m;
  for (m = 0; m < expo->dim; ++m) pin[m] = 1.0*rand()/RAND_MAX;

  sir_init_index(expo);
  sir_init_all(expo);
  sir_init_lambda(expo);
  init_mupsi_const(expo);

  /* setup full matrix */
  int a, b;
  double* A = (double*) calloc(expo->dim*expo->dim,sizeof(double));
  for (m = 0; m < expo->dim; ++m) {
    I = expo->lookup[m];
    R = expo->lookup[m+expo->dim];
    a = sir_index(I-1,R+1,expo);
    b = sir_index(I+1, R ,expo);
    A[m*expo->dim+m] = expo->mat[m];
    if (a >= 0) A[a*expo->dim+m] = expo->mat[m+expo->dim];
    if (b >= 0) A[b*expo->dim+m] = expo->mat[m+2*expo->dim];
  }

  for (m = 0; m < expo->dim; ++m) {
    int l;
    for (l = 0; l < expo->dim; ++l) {
      if (A[l*expo->dim+m] != 0.0) printf("#");
      else printf(".");
    }
    printf("\n");
  }

  double nrm = dnrm2_(&expo->dim,pin,&one);
  printf("nrm = %f\n",nrm);

  char trans = 't';
  double alpha = 1.0;
  double dzero = 0.0;
  sir_matfunc(trans,0,0,alpha,pin,pout,expo);

  double nrm_out = dnrm2_(&expo->dim,pout,&one);
  printf("nrm = %f\n",nrm_out);

  dgemv_(&trans,&expo->dim,&expo->dim,&alpha,A,&expo->dim,pin,&one,&dzero,pout2,&one);

  nrm_out = dnrm2_(&expo->dim,pout,&one);
  printf("nrm(full) = %f\n",nrm_out);

  printf("%d\n",expo->lookup[expo->dim+3]);

  m = 0;
  for (m = 0; m < expo->dim; ++m) {
    I = expo->lookup[m];
    R = expo->lookup[m+expo->dim];
    printf("(%2d|%2d,%2d) | % 6f | % 6f | % 6f\n",m,I,R,pin[m],pout[m],pout2[m]);
  }

  free(pin);
  free(pout);
  free(pout2);
  free(A);

  expo_free_all(expo);

  return 0;
}


