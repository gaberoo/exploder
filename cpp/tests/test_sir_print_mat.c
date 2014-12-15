#include <stdlib.h>

#include "shared/sir.h"

int main() {

  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));
  init_expo(expo);

  sir_start(expo,3);
  printf("dim = %d\n",expo->dim);

  expo_alloc_all(expo);
  sir_init_all(expo);
  sir_init_index(expo);

  int m;
  int a, b;
  int I, R;

  for (m = 0; m < expo->dim; ++m) {
    I = expo->lookup[m];
    R = expo->lookup[m+expo->dim];
    printf("%d --> (%d,%d)\n",m,I,R);
  }

  printf("-------\n");

  for (m = 0; m < expo->dim; ++m) {
    I = expo->lookup[m];
    R = expo->lookup[m+expo->dim];
    a = sir_index(I+1,R,expo);
    b = sir_index(I-1,R+1,expo);
    printf("p'(%d) =",m);
    printf(" A(%d,%d|%d:%d) p(%d)",m,m,I,R,m);
    if (I+R < expo->N_max) printf(" + A(%d,%d|%d:%d) p(%d)",m,a,I+1,R,a);
    if (I > 0) printf(" + A(%d,%d|%d:%d) p(%d)",m,b,I-1,R+1,b);
    printf("\n");
  }

  printf("-------\n");

  for (m = 0; m < expo->dim; ++m) {
    I = expo->lookup[m];
    R = expo->lookup[m+expo->dim];
    a = sir_index(I-1,R,expo);
    b = sir_index(I+1,R-1,expo);
    printf("p'(%d) =",m);
    printf(" A(%d,%d|%d:%d) p(%d)",m,m,I,R,m);
    if (R > 0) printf(" + A(%d|%d:%d,%d) p(%d)",b,I+1,R-1,m,b);
    if (I > 0) printf(" + A(%d|%d:%d,%d) p(%d)",a,I+1,R-1,m,a);
    printf("\n");
  }

  expo_free_all(expo);

  return 0;
}

