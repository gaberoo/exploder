#include "shared/sir.h"

int main(int argc, char** argv) {
  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);

  expo->N_max = (argc > 1) ? atoi(argv[1]) : 10;
  expo->dim = sir_index(0,expo->N_max,expo);
  expo->lookup = (int*) malloc(2*expo->dim*sizeof(int));
  sir_init_index(expo);

  int I, R;
  int m = 0;
  for (R = 0; R <= expo->N_max; ++R) {
    for (I = 0; I <= expo->N_max-R; ++I) {
      printf("%2d --> (%2d,%2d) ",m++,I,R);
    }
    printf("\n");
  }

  for (R = 0; R <= expo->N_max; ++R) {
    for (I = 0; I <= expo->N_max-R; ++I) {
      printf("%4d    ",sir_index(I,R,expo));
    }
    printf("\n");
  }

  for (m = 0; m < expo->dim; ++m) {
    I = expo->lookup[m];
    R = expo->lookup[m+expo->dim];
    printf("%2d --> (%2d,%2d)\n",m,I,R);
  }

  printf("Max index = %d\n",sir_index_N(0,expo->N_max,expo->N_max));

  free(expo->lookup);
  expo_type_free(expo);

  return 0;
}

