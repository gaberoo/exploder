#pragma OPENCL EXTENSION cl_khr_fp64 : enable

/* from: AMD -- OpenCL(TM) Optimization Case Study: Simple Reductions */
__kernel void inf_norm_kernel(__global clReal* buffer,
                              __local clReal* scratch,
                              const int len,
                              __global clReal* result) 
{
  int global_index = get_global_id(0);
  int local_index = get_local_id(0);

  clReal accumulator = 0.0;
  clReal other;
  clReal mine;
  clReal element;

  // Loop sequentially over chunks of input vector
  while (global_index < len) {
    element = fabs(buffer[global_index]);
    if (accumulator < element) accumulator = element;
    global_index += get_global_size(0);
  }

  // Perform parallel reduction
  scratch[local_index] = accumulator;
  barrier(CLK_LOCAL_MEM_FENCE);
  for (int offset = get_local_size(0) / 2; offset > 0; offset >>= 1) {
    if (local_index < offset) {
      other = scratch[local_index + offset];
      mine  = scratch[local_index];
      if (other > mine) scratch[local_index] = other;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (local_index == 0) result[get_group_id(0)] = scratch[0];
}

/****************************************************************************/

__kernel void two_norm_kernel(__global clReal* buffer,
                              __local clReal* scratch,
                              const int len,
                              __global clReal* result) 
{
  int global_index = get_global_id(0);
  int local_index = get_local_id(0);

  /* Loop sequentially over chunks of input vector */
  clReal accumulator = 0.0f;
  while (global_index < len) {
    accumulator += buffer[global_index]*buffer[global_index];
    global_index += get_global_size(0);
  }
  scratch[local_index] = accumulator;
  barrier(CLK_LOCAL_MEM_FENCE);

  /* Perform parallel reduction */
  for (int offset = get_local_size(0)/2; offset > 0; offset >>= 1) {
    if (local_index < offset) {
      scratch[local_index] += scratch[local_index + offset];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (local_index == 0) result[get_group_id(0)] = scratch[0];
//  if (local_index == 0) {
//    result[0] += scratch[0];
//    // barrier(CLK_GLOBAL_MEM_FENCE);
//  }
}

/****************************************************************************/

__kernel void one_norm_tridiag_mat_kernel(__global clReal* mat,
                                          __local clReal* scratch,
                                          const int len,
                                          const int lda,
                                          __global clReal* result) 
{
  int global_index = get_global_id(0);

  /* Loop sequentially over chunks of input vector */
  clReal accumulator = 0.0f;

  clReal element;
  while (global_index < len) {
    element  = fabs(mat[global_index]);
    element += fabs(mat[lda+global_index]);
    element += fabs(mat[2*lda+global_index]);

    accumulator = (accumulator > element) ? accumulator : element;

    global_index += get_global_size(0);
  }

  /* Perform parallel reduction */
  int local_index = get_local_id(0);

  scratch[local_index] = accumulator;
  barrier(CLK_LOCAL_MEM_FENCE);

  for (int offset = get_local_size(0)/2; offset > 0; offset >>= 1) {
    if (local_index < offset) {
      clReal other = fabs(scratch[local_index + offset]);
      clReal mine  = fabs(scratch[local_index]);
      scratch[local_index] = (mine > other) ? mine : other;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (local_index == 0) result[get_global_id(0)] = scratch[0];
}

/****************************************************************************/

__kernel void tridiag_kernel (__global const clReal* A,
                              __global const clReal* x,
                              __global clReal* res,
                              const char trans,
                              const clReal alpha, 
                              const int N)
{
  const int row = get_global_id(0);
  if (row < N) {
    clReal out;
    switch (trans) {
      case 'n':
      case 'N':
      default:
        out = A[N+row] * x[row];
        if (row > 0) out += A[row] * x[row-1];
        if (row < N-1) out += A[2*N+row] * x[row+1];
        break;
      case 't':
      case 'T':
        out = A[N+row] * x[row];
        if (row > 0) out += A[2*N+row-1] * x[row-1];
        if (row < N-1) out += A[row+1] * x[row+1];
        break;
    }
    res[row] = alpha*out;
  }
}

/****************************************************************************/

__kernel void vector_add_kernel (__global clReal* a,
                                 __global const clReal* b,
                                 const int N)
{
  const int row = get_global_id(0);
  if (row < N) a[row] += b[row];
}

/****************************************************************************/

__kernel void vector_copy_kernel (__global clReal* a,
                                  __global const clReal* b,
                                  const int N)
{
  const int row = get_global_id(0);
  if (row < N) a[row] = b[row];
}

/****************************************************************************/

__kernel void vector_smult_kernel (__global clReal* x,
                                   const clReal c,
                                   const int N)
{
  const int row = get_global_id(0);
  if (row < N) x[row] *= c;
}

/****************************************************************************/

__kernel void expo_samp_kernel(__global const clReal* pin,
                               __global clReal* pout,
                               __global const clReal* psi,
                               const int N)
{
  const int m = get_global_id(0);
  if (m < N) {
    if (m == 0) pout[m] = 0.0f;
    else pout[m] = psi[m]*pin[m-1];
  }
}

/****************************************************************************/

__kernel void expo_trans_kernel(__global const clReal* pin,
                                __global clReal* pout,
                                __global const clReal* lambda,
                                const int N)
{
  const int m = get_global_id(0);
  if (m < N-1) pout[m] = 2.0f*lambda[m]*pin[m+1];
  else if (m == N-1) pout[m] = 0.0f;
}

/****************************************************************************/

__kernel void init_mat_kernel(__global clReal* lambda,
                              __global clReal* mat,
                              __global clReal* mu,
                              __global clReal* psi,
                              const int ki,
                              const clReal shift,
                              const int N) 
{
  const int m = get_global_id(0);
  if (m < ki) {
    mat[m] = 0.0f;
    mat[m+N] = 0.0f;
    mat[m+2*N] = 0.0f;
  } else if (m < N) {
    mat[m]     = (m > ki) ? (m-ki)*mu[m] : 0.0f;      /* lower diagonal */
    mat[N+m]   = -m*(lambda[m]+psi[m]+mu[m]) - shift; /*    diagonal    */
    mat[2*N+m] = (m < N-1) ? (m+ki)*lambda[m] : 0.0f; /* upper diagnoal */
  }
}

/****************************************************************************/

__kernel void neg_vals_kernel(__global clReal* buffer,
                              __global clReal* out,
                              __local int* scratch,
                              const int len,
                              __global clReal* result) 
{
  int global_index = get_global_id(0);
  int local_index = get_local_id(0);

  /* Loop sequentially over chunks of input vector */
  int accumulator = 0;
  while (global_index < len) {
    clReal element = buffer[global_index];
    if (element < 0.0f) {
      out[global_index] = 0.0f;
      ++accumulator;
    } else {
      out[global_index] = element;
    }
    global_index += get_global_size(0);
  }
  scratch[local_index] = accumulator;
  barrier(CLK_LOCAL_MEM_FENCE);

  /* Perform parallel reduction */
  for (int offset = get_local_size(0)/2; offset > 0; offset >>= 1) {
    if (local_index < offset) {
      scratch[local_index] += scratch[local_index + offset];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (local_index == 0) result[get_group_id(0)] = scratch[0];
}

/****************************************************************************/
/**************************** INDEX KERNELS *********************************/
/****************************************************************************/

int SIR(int I, int R, int N) {
  if (I >= 0 && R >= 0 && I+R <= N && I <= N && R <= N) {
    return R*(N+1) - R*(R-1)/2 + I;
  } else {
    return -1;
  }
}

__kernel void tridiag_sir_kernel (__global const clReal* A,
                                  __global const clReal* x,
                                  __global const int* index,
                                  __global clReal* res,
                                  const char trans,
                                  const clReal alpha, 
                                  const int N,
                                  const int dim)
{
  const int m = get_global_id(0);
  if (m < dim) {
    int I = index[m];               /* the 'index' array stores the ... */
    int R = index[m+dim];           /* ... (I,R) values */
    int a, b;

    clReal out;
    out = A[m]*x[m];

    switch (trans) {
      case 'n':
      case 'N':
      default:
        a = SIR(I-1,R+1,N);
        b = SIR(I+1,R,N);
        if (a >= 0) out += A[  dim+m] * x[a];
        if (b >= 0) out += A[2*dim+m] * x[b];
        break;
      case 't':
      case 'T':
        a = SIR(I+1,R-1,N);
        b = SIR(I-1,R,N);
        if (a >= 0) out += A[  dim+a] * x[a];
        if (b >= 0) out += A[2*dim+b] * x[b];
        break;
    }

    res[m] = alpha*out;
  }
}

/****************************************************************************/

__kernel void expo_sir_samp_kernel(__global const clReal* pin,
                                   __global clReal* pout,
                                   __global const int* lookup,
                                   __global const clReal* psi,
                                   const int N,
                                   const int dim)
{
  const int m = get_global_id(0);
  if (m < dim) {
    int I = lookup[m];
    int R = lookup[m+dim];
    int a = SIR(I-1,R+1,N);
    if (a >= 0 && I > 0) pout[m] = psi[m]*pin[a];
    else pout[m] = 0.0f;
  }
}

/****************************************************************************/

__kernel void expo_sir_trans_kernel(__global const clReal* pin,
                                    __global clReal* pout,
                                    __global const int* lookup,
                                    __global const clReal* lambda,
                                    const int N,
                                    const int dim)
{
  const int m = get_global_id(0);
  if (m < dim) {
    int I = lookup[m];
    int R = lookup[m+dim];
    int a = SIR(I+1,R,N);
    if (a >= 0 && I < N-R) pout[m] = 2.0f*lambda[m]*pin[a];
    else pout[m] = 0.0f;
  }
}

/****************************************************************************/


