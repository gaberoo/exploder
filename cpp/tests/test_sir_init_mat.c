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

  int dim = sir_index(0,expo->N_max,expo)+1;
  printf("dim = %d\n",dim);

  expo->lambdaVec = (double*) malloc(dim*sizeof(double));
  expo->mat = (double*) malloc(3*dim*sizeof(double));
  expo->mat_i = (int*) malloc(3*dim*sizeof(int));

  sir_init_all(expo);

  int m = 0;
  for (m = 0; m < dim; ++m) {
    printf("%2d:% 6.2f | %2d:% 6.2f | %2d:% 6.2f\n",
           expo->mat_i[m],expo->mat[m],
           expo->mat_i[dim+m],expo->mat[dim+m],
           expo->mat_i[2*dim+m],expo->mat[2*dim+m]);
  }

  free(expo->lambdaVec);
  free(expo->mat_i);
  free(expo->mat);

  return 0;
}


