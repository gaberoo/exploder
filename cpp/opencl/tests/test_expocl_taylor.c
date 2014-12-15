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

  expo_cl* ecl = expocl_init(expo,CL_DEVICE_TYPE_GPU,128);
  expocl_init_kernel(ecl,expo);
  expocl_create_buffers(ecl);
  expocl_copy_mat(ecl,expo);

  int n = N;
  clReal* cl_b1 = (clReal*) malloc(n*sizeof(clReal));
  for (i = 0; i < n; ++i) cl_b1[i] = (clReal) x[i];

  expocl_copy_x(ecl,cl_b1);

  /* copy b1 onto b on device */
  expocl_copy_bx(ecl);

  int s = 10;
  int tcol = 4;
  double t = 1.0;
  int mv = 0;
  double c1, c2;
  int full_term = 0;
  double bnorm;
  double tol = 1e-5;
  clReal eta = 0.2;

  printf("Go!\n");

  int k;
  int ss;
  for (ss = 1; ss <= s; ++ss) {
    /* get norm of input vector/matrix */
    c1 = expocl_inf_norm_x(ecl);

    for (k = 1; k <= tcol; ++k) {
      /* calculate y = A x */
      expocl_calc(ecl,t/(s*k)); 

      /* x => y */
      expocl_swap_xy(ecl);
      ++mv;  /* increment matrix-vector product counter */

      /* norm of output vector */
      c2 = expocl_inf_norm_x(ecl);

      /* add x onto b */
      expocl_add_bx(ecl);

      /* check whether convergence has been reached before full truncation */
      if (! full_term) {
        /* get new infinite norm of the new b */
        bnorm = expocl_inf_norm_b(ecl);
        
        if (c1+c2 <= tol*bnorm) break;

        c1 = c2;
      }
      
      printf("k = %d: c1 = %g, c2 = %g\n",k,c1,c2);
    }

    expocl_smult_b(ecl,eta);
    expocl_copy_xb(ecl);
  }

  free(cl_b1);
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

