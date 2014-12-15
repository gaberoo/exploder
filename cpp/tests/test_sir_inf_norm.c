#include <stdlib.h>

#include "shared/sir.h"

int main() {
  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);

  expo->N_max = 1000;
  expo->K = 1000.0;
  expo->beta = 1.0;
  expo->mu = 0.2;
  expo->psi = 0.2;

  int dim = sir_index(0,expo->N_max,expo)+1;
  printf("dim = %d\n",dim);

  expo->mat = (double*) malloc(3*dim*sizeof(double));
  expo->mat_i = (int*) malloc(3*dim*sizeof(int));

  sir_init_mat(expo);

  double inf_norm = sir_inf_norm(expo);

  printf("|A|_inf = %g\n",inf_norm);

  free(expo->mat_i);
  free(expo->mat);

  return 0;
}


