#include <stdlib.h>
#include <stdio.h>

double dlange_(char* norm, int* M, int* N, double* A, int* LDA, void* WORK);
 
int main() {

  int M = 3;
  int N = 2;
  double A[] = { 1., 4., -3., -10., 6., 4.5 };
  void* work = malloc(M*sizeof(double));

  char norm;

  norm = 'o';
  double oneNorm = dlange_(&norm,&M,&N,A,&M,work);

  norm = 'i';
  double infNorm = dlange_(&norm,&M,&N,A,&M,work);

  int i, j;
  for (i = 0; i < M; ++i) {
    printf("| ");
    for (j = 0; j < N; ++j) {
      printf("%8.2f ",A[j*M+i]);
    }
    printf(" |\n");
  }
  printf("one-norm = %g\n",oneNorm);
  printf("inf-norm = %g\n",infNorm);

  free(work);
  return 0;
}

