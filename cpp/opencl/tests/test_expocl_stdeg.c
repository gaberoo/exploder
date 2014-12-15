#include <stdlib.h>
#include <stdio.h>

#include "../shared/expo.h"
#include "../shared/expmv.h"

#include "expocl.h"
#include "cl_expomv.h"

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

  expo->init_all(expo);

  int i = 0;

  /********************** expocl ************************/

  int info = 0;

  /* calculate required memory and allocate */
  int ncol = 1;
  int memlen = (2*expo->dim+(expo->p_max-1)*expo->m_max);
  int wlen   = 2*expo->p_max + (6*ncol+3)*expo->dim + ncol + expo->m_max*(expo->p_max-1);
  int iwlen  = 2*expo->dim + 4;

  double* wrk = (double*) malloc((memlen+wlen)*sizeof(double));
  int* iwrk = (int*) malloc(iwlen*sizeof(int));

  double* p0 = wrk;
  double* pT = p0 + N;
  double* tm = pT + N;

  /* set present time */
  wrk[0] = 0.0;

  clReal* b = (clReal*) malloc(expo->dim*sizeof(clReal));
  for (i = 0; i < expo->dim; ++i) {
    b[i] = 1.0;
    p0[i] = b[i];
  }

  double dt = 1.0;

  print_info();

  expo_cl* ecl = expocl_init(expo,CL_DEVICE_TYPE_GPU,128);
  // expo_cl* ecl = expocl_init(expo,CL_DEVICE_TYPE_CPU,1);
  expocl_init_kernel(ecl,expo);
  expocl_create_buffers(ecl);
  expocl_copy_mat(ecl,expo);

  expocl_copy_x(ecl,b);
  expocl_copy_bx(ecl);

  double mdt = max_dt(expo);
  printf("max dt = %g\n",mdt);

  int est_norm = 0;
  double* talpha = (double*) calloc(expo->p_max,sizeof(double));
  double* teta   = (double*) calloc(expo->p_max,sizeof(double));

  dt = 0.5*mdt;
  int mv = cl_stdeg(dt, b, expo, ecl, est_norm, talpha, teta, tm,
                    0, NULL, 0, NULL);
  printf("mv = %d\n",mv);
  for (i = 0; i < expo->p_max; ++i) {
    printf("alpha = %g, eta = %g\n",talpha[i],teta[i]);
  }

  dt = 1.5*mdt;
  mv = cl_stdeg(dt, b, expo, ecl, est_norm, talpha, teta, tm,
                    0, NULL, 0, NULL);
  printf("mv = %d\n",mv);

  free(teta);
  free(talpha);

  expocl_free_kernel(ecl);
  expocl_free(ecl);

  printf(" info = %d\n",info);

  free(wrk);
  free(iwrk);

  /********************** clean up **********************/

  free(expo->mat);
  expo->mat = NULL;

  free(expo->lambdaVec);
  expo->lambdaVec = NULL;

  expo_type_free(expo);

  free(b);

  return 0;
}

