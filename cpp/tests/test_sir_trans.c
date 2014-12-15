#include <stdlib.h>

#include "shared/sir.h"

int main() {
  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);

  expo->N_max = 10;
  expo->K = 10.0;
  expo->beta = 1.0;
  expo->mu = 0.2;
  expo->psi = 0.2;
  expo->ki = 3;

  int dim = sir_index(0,expo->N_max,expo)+1;
  printf("dim = %d\n",dim);

  expo->lambdaVec = (double*) malloc(dim*sizeof(double));
  expo->mat = (double*) malloc(3*dim*sizeof(double));
  expo->mat_i = (int*) malloc(3*dim*sizeof(int));

  double* pin  = (double*) calloc(dim,sizeof(double));
  double* pout = (double*) calloc(dim,sizeof(double));
  double* pout2 = (double*) calloc(dim,sizeof(double));

  int m = 0;
  for (m = 0; m < dim; ++m) pin[m] = 1.0*m/dim;

  sir_init_all(expo);
  sir_trans(pin,pout,expo);
  sir_sample(pin,pout2,expo);

  int I, R;
  for (R = 0; R <= expo->N_max; ++R) {
    for (I = 0; I <= expo->N_max-R; ++I) {
      m = sir_index(I,R,expo);
      printf("%2d | %2d,%2d | %8g | %6.2f %6.2f %6.2f\n",
             m,I,R,expo->lambdaVec[m],pin[m],pout[m],pout2[m]);
    }
  }

  free(pin);
  free(pout);
  free(pout2);

  free(expo->lambdaVec);
  free(expo->mat_i);
  free(expo->mat);

  return 0;
}


