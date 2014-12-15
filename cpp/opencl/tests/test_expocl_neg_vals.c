#include <stdlib.h>
#include <stdio.h>

#include "../shared/expo.h"
#include "../shared/expmv.h"

#include "expocl.h"

double max_dt(expo_type* expo);
int cl_expmv(double t, clReal* b, double* tm, int recalcm,
             char prec, int full_term, int prnt,
             int est_norm, expo_type* expo, expo_cl* ecl,
             int wrklen, double* wrk, int iwrklen, int* iwrk);

float snrm2_(int* n, float* x, int* incx);

int main(int argc, char** argv) {
  double RN[]    = { 100 };
  double Rbeta[] = { 1.0 };
  double Rmu[]   = { 20.0 };
  double Rpsi[]  = { 0.1 };
  double Rki[] = { 10 };
  int estimateNorm[] = { 0 };

  if (argc > 2) RN[0] = atoi(argv[2]);

  Rbeta[0] /= RN[0];
  Rmu[0] /= RN[0];
  Rpsi[0] /= RN[0];

  expo_type* expo = expo_type_alloc();

  expo->N       = (int) ceil(RN[0]); /* total population size */
  expo->K       = RN[0];             /* carrying capacity */
  expo->ki      = *Rki;              /* current number of lineages */
  expo->beta    = Rbeta[0];          /* branching rate */
  expo->mu      = Rmu[0];            /* extinction rate */
  expo->psi     = Rpsi[0];           /* sampling rate */
  expo->SImodel = 1;                 /* Model: DD (SI=1), INF (SI=0) */
  expo->rescale = 1;                 /* rescale probability vector */
  expo->vflag   = 3;                 /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = 2;               /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = RN;                /* N parameters */
  expo->betaVec = Rbeta;             /* beta parameters */
  expo->muVec   = Rmu;               /* mu parameters */
  expo->psiVec  = Rpsi;              /* psi parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */
  expo->est_norm = *estimateNorm;    /* force estimation of matrix norm */

  expo->norm = &matFuncOneNormAnalytical;

  expo->N_max = RN[0];
  int N = RN[0] + 1;
  expo->dim = N;

  /* allocate workspace for matrix */
  expo->mat = (double*) malloc(3*N*sizeof(double));
  expo->lambda = &lambdaSI;

  /* allocate workspace for lambda */
  expo->lambdaVec = (double*) malloc(N*sizeof(double));
  expo->muFun = (double*) malloc(N*sizeof(double));
  expo->psiFun = (double*) malloc(N*sizeof(double));

  int i = 0;
  for (i = 0; i < N; ++i) {
    expo->muFun[i] = expo->mu;
    expo->psiFun[i] = expo->psi;
  }

  expo->init_all(expo);

  /********************** expocl ************************/

  double* p0 = (double*) malloc(expo->dim*sizeof(double));
  clReal* b = (clReal*) malloc(expo->dim*sizeof(clReal));

  for (i = 0; i < expo->dim; ++i) {
    b[i] = -1.0;
    p0[i] = b[i];
  }

  int one = 1;

  print_info();

  printf("Initializing device...");
  expo_cl* ecl = expocl_init(expo,CL_DEVICE_TYPE_GPU,128);
  printf("done.\n");

  printf("Initializing kernels and buffers...");
  expocl_init_kernel(ecl,expo);
  expocl_create_buffers(ecl);
  expocl_create_vector_buffers(ecl);
  printf("done.\n");

  expocl_copy_mat(ecl,expo);
  expocl_copy_vectors(ecl,expo);

  expocl_copy_x(ecl,b);
  expocl_copy_bx(ecl);

#ifndef DOUBLE
  expocl_read_b(ecl,b);
  float s_norm = snrm2_(&N,b,&one);
  printf("|p|_s = %g\n",s_norm);
#endif

  clReal gpu_norm = expocl_two_norm_b(ecl);
  double norm = dnrm2_(&N,p0,&one);
  printf("|p|   = %g\n",norm);
  printf("|p|_g = %g\n",gpu_norm);

  int found_neg = expocl_neg_vals(ecl);
  printf("# negative = %d\n",found_neg);

  expocl_read_x(ecl,b);
  for (i = 0; i < 20; ++i) printf("%d %g\n",i,b[i]);

  expocl_release_vector_buffers(ecl);
  expocl_release_buffers(ecl);
  expocl_free_kernel(ecl);
  expocl_free(ecl);

  /********************** clean up **********************/

  free(expo->mat);
  free(expo->lambdaVec);
  free(expo->muFun);
  free(expo->psiFun);

  expo_type_free(expo);

  free(b);
  free(p0);

  return 0;
}

