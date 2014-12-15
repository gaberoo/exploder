#include <stdlib.h>
#include <stdio.h>

#include "../shared/expo.h"
#include "../shared/expmv.h"

#include "expocl.h"
#include "cl_expomv.h"

#include "tests/start_vec.h"

float snrm2_(int* n, float* x, int* incx);

int main(int argc, char** argv) {
  double RN[]    = { 30 };
  double Rbeta[] = { 1.0 };
  double Rmu[]   = { 0.1 };
  double Rpsi[]  = { 0.1 };
  double Rki[] = { 1 };
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
  expo->SImodel = 2;                 /* Model: DD (SI=1), INF (SI=0) */
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

  switch (expo->SImodel) {
    case 1:
    case 0:
    default:
      /* SIS */
      expo_start(expo,max_pop_size(expo));
      expo->init_all = &init_all;
      expo->matvec = &matFuncExpmv;
      break;

    case 2:
      /* SIR */
      sir_start(expo,max_pop_size(expo));
      expo->init_all = &sir_init_all;
      expo->matvec   = &sir_matfunc;
      expo->ft       = &sir_trans;
      expo->fs       = &sir_sample;
      expo->norm     = &sir_one_norm;
      expo->trace    = &sir_trace;
      break;
  }

  expo_alloc_all(expo);
  expo->init_all(expo);

  sir_init_lambda(expo);
  init_mupsi_const(expo);
  init_wrk(expo);

  /********************** expocl ************************/

  /* calculate required memory and allocate */

  double* p0 = expo->wrk;
  double* pT = p0 + expo->dim;
  double* start_vec = pT + expo->dim;
  double* wrk = start_vec + expo->dim;

  // for (i = 0; i < expo->dim; ++i) p0[i] = 1.0;
  memcpy(p0,start_vec,expo->dim*sizeof(double));

  double* b = (double*) malloc(expo->dim*sizeof(double));
  memcpy(b,p0,expo->dim*sizeof(double));

  expo_cl* ecl = expocl_init(expo,CL_DEVICE_TYPE_GPU,128);

  expocl_init_kernel(ecl,expo);

  expocl_create_buffers(ecl);
  expocl_create_index_buffer(ecl,expo);

  expocl_copy_mat(ecl,expo);
  expocl_copy_index(ecl,expo);

  printf("\n");

  double est = 0.0;
  int mv = 0;
  expocl_normam(1.0,1,&est,&mv,wrk,expo->iwrk,expo,ecl);
  printf("|A| = %f\n",est);
  printf("mv  = %d\n",mv);

  mv = normAm(1.0,1,&est,wrk,expo->iwrk,expo);
  printf("|A| = %f\n",est);
  printf("mv  = %d\n",mv);

  expocl_free_kernel(ecl);
  expocl_free(ecl);

  free(b);

  /********************** clean up **********************/

  expo_free_all(expo);
  expo_type_free(expo);

  return 0;
}

