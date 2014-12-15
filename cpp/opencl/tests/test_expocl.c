#include <stdlib.h>
#include <stdio.h>

#include "expocl.h"

void tridiag_seq (const double* A, const double* x, double* res, 
                  clReal alpha, const int N)
{
  int row = 0;

  clReal cc = 0.0;
  clReal aa =   A[N+row]*x[row];
  clReal bb = (N > 1) ? A[2*N+row]*x[row+1] : 0.0;
  res[row] = alpha*(aa+bb+cc);

  for (row = 1; row < N-1; ++row) {
    cc =     A[row]*x[row-1];
    aa =   A[N+row]*x[row];
    bb = A[2*N+row]*x[row+1];
    res[row] = alpha*(aa+bb+cc);
  }

  if (N > 1) {
    row = N-1;
    cc =     A[row]*x[row-1];
    aa =   A[N+row]*x[row];
    bb = 0.0;
    res[row] = alpha*(aa+bb+cc);
  }
}

int main(int argc, char** argv) {
  double RN[]    = { 1000 };
  double Rbeta[] = { 1.0 };
  double Rmu[]   = { 0.2 };
  double Rpsi[]  = { 0.001 };
  double Rki[] = { 10 };
  int estimateNorm[] = { 1 };

  expo_type* expo = expo_type_alloc();

  expo->N       = (int) ceil(RN[0]); /* total population size */
  expo->K       = RN[0];             /* carrying capacity */
  expo->ki      = *Rki;              /* current number of lineages */
  expo->beta    = Rbeta[0];          /* branching rate */
  expo->mu      = Rmu[0];            /* extinction rate */
  expo->psi     = Rpsi[0];           /* sampling rate */
  expo->SImodel = 1;                 /* Model: DD (SI=1), INF (SI=0) */
  expo->rescale = 1;                 /* rescale probability vector */
  expo->vflag   = 0;                 /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = 2;               /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = RN;                /* N parameters */
  expo->betaVec = Rbeta;             /* beta parameters */
  expo->muVec   = Rmu;               /* mu parameters */
  expo->psiVec  = Rpsi;              /* psi parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */
  expo->est_norm = *estimateNorm;    /* force estimation of matrix norm */

  expo->N_max = RN[0];
  int N = RN[0] + 1;

  /* allocate workspace for matrix */
  expo->mat = (double*) malloc(3*N*sizeof(double));
  expo->lambda = &lambdaSI;

  /* allocate workspace for lambda */
  expo->lambdaVec = (double*) malloc(N*sizeof(double));
  init_lambda(expo);

  init_mat(expo);

  int i = 0;
//  for (i = 0; i < N; ++i) {
//        expo->mat[i] = 0.0;
//      expo->mat[N+i] = 1.0;
//    expo->mat[2*N+i] = 0.0;
//  }

  /********************** expocl ************************/

  clReal* x = (clReal*) malloc(N*sizeof(clReal));
  clReal* y = (clReal*) calloc(N,sizeof(clReal));

  for (i = 0; i < N; ++i) x[i] = 1.0; 

  double* dx = (double*) calloc(N,sizeof(double));
  double* dy = (double*) calloc(N,sizeof(double));
  for (i = 0; i < N; ++i) dx[i] = (double) x[i];

  clReal alpha = 0.3;
  clReal norm = 0.0;

  expo_cl* ecl = expocl_init(expo,CL_DEVICE_TYPE_GPU,128);
  expocl_init_kernel(ecl,expo);
  expocl_create_buffers(ecl);
  expocl_init_mat(ecl,expo);

  expocl_copy_x(ecl,x);
  expocl_copy_bx(ecl);
  expocl_read_b(ecl,x);
  printf("b(0) = %g\n",x[0]);

  norm = expocl_inf_norm_b(ecl);
  printf("norm = %g\n",norm);

  expocl_calc(ecl,alpha); expocl_swap_xy(ecl);
  expocl_calc(ecl,alpha); expocl_swap_xy(ecl);
  expocl_calc(ecl,alpha); expocl_swap_xy(ecl);
  
  expocl_read_x(ecl,y);

  norm = expocl_inf_norm_x(ecl);
  printf("norm = %g\n",norm);

  expocl_release_buffers(ecl);
  expocl_free_mat(ecl);

  tridiag_seq(expo->mat,dx,dy,alpha,N);
  tridiag_seq(expo->mat,dy,dx,alpha,N);
  tridiag_seq(expo->mat,dx,dy,alpha,N);

  // for (i = 0; i < N; ++i) 
  for (i = N/2-10; i < N/2+10; ++i) {
    printf("%4d | %12g | %12g \n",i,y[i],dy[i]);
  }

  expocl_free(ecl);

  free(dx);
  free(dy);

  free(x);
  free(y);

  /********************** clean up **********************/

  free(expo->mat);
  expo->mat = NULL;

  free(expo->lambdaVec);
  expo->lambdaVec = NULL;

  free(expo);
  expo = NULL;

  return 0;
}

