#include <stdlib.h>

#include "shared/sir.h"

int main() {
  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);

  expo->N_max = 4;
  expo->K = 4.0;
  expo->beta = 1.0;
  expo->mu = 0.2;
  expo->psi = 0.2;

  double S;
  double lambda;

  int dim = sir_index(0,expo->N_max,expo);

  double* pin  = (double*) calloc(dim,sizeof(double));
  double* pout = (double*) calloc(dim,sizeof(double));

  int i, j1, j2;
  /* This is acutally a linear loop ! */
  int I, R;
  for (R = 0; R <= expo->N_max; ++R) {
    for (I = 0; I <= expo->N_max-R; ++I) {
      S = expo->K-R-I;

      if (S < 0.0) S = 0.0;
      lambda = expo->beta*I*S;

      double a = -(lambda-expo->mu-expo->psi);
      i = sir_index(I,R,expo);
      pout[i] = a*pin[i];

      if (I > 0 && R < expo->N_max) {
        j1 = sir_index(I-1,R+1,expo);
        pout[i] += (I-expo->ki)*expo->mu*pin[j1];
      } else {
        j1 = -1;
      }

      if (I+R < expo->N_max) {
        j2 = sir_index(I+1,R,expo);
        pout[i] += (I+expo->ki)*lambda*pin[j2];
      } else {
        j2 = -1;
      }

      printf("%2d (%2d,%2d) | %2d (%2d,%2d) | %2d (%2d,%2d)\n",
             i,I,R,j1,I-1,R+1,j2,I+1,R);

      // printf("%2d | %2d | %2d \n",i,j1,j2);

    }
  }

  free(pin);
  free(pout);

  return 0;
}


