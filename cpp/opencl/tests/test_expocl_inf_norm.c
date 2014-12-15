#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int main() {
  int N = 11;
  double* vec = (double*) malloc(N*sizeof(double));

  int i;
  for (i = 0; i < N; ++i) {
    vec[i] = (i-N/2)*(i-N/2);
    printf("%f ",vec[i]);
  }
  printf("\n");

  int global_size = 2;
  int local_size = N/global_size;
  printf("local_size = %d\n");

  int local_index;
  int wrk_item;

  for (wrk_item = 0; wrk_item < global_size; ++wrk_item) {
    printf("WRK GRP = %d\n",wrk_item);
    for (local_index = 0; local_index < local_size; ++local_index) {
      printf("  LOCAL = %d\n",local_index);
      int gi = wrk_item*global_size + local_index;
      double accumulator = 0.0f;
      double other;
      double mine;
      double element;
      // Loop sequentially over chunks of input vector
      while (gi < N) {
        element = fabs(vec[gi]);
        printf("    %d %f\n",gi,element);
        if (element > accumulator) accumulator = element;
        gi += global_size;
      }
    }
  }

//  // Perform parallel reduction
//  scratch[local_index] = accumulator;
//  barrier(CLK_LOCAL_MEM_FENCE);
//  for (int offset = get_local_size(0) / 2; offset > 0; offset >>= 1) {
//    if (local_index < offset) {
//      other = scratch[local_index + offset];
//      mine  = scratch[local_index];
//      if (other > mine) scratch[local_index] = other;
//    }
//    barrier(CLK_LOCAL_MEM_FENCE);
//  }
//
//  if (local_index == 0) result[get_global_id(0)] = scratch[0];

  free(vec);

  return 0;
}

