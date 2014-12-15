#include <stdlib.h>
#include <stdio.h>

#include "../shared/expo.h"
#include "../shared/expmv.h"

#include "expocl.h"
#include "cl_expomv.h"

#include <time.h>
#include <mach/clock.h>
#include <mach/mach.h>

float snrm2_(int* n, float* x, int* incx);

int main(int argc, char** argv) {
  double RN[]    = { 2000 };
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
  init_expo(expo);

  expo->SImodel = (argc > 1) ? atoi(argv[1]) : 1; /* Model: DD (SI=1), INF (SI=0) */
  expo->N       = (int) ceil(RN[0]); /* total population size */
  expo->K       = RN[0];             /* carrying capacity */
  expo->ki      = *Rki;              /* current number of lineages */
  expo->beta    = Rbeta[0];          /* branching rate */
  expo->mu      = Rmu[0];            /* extinction rate */
  expo->psi     = Rpsi[0];           /* sampling rate */
  expo->rescale = 1;                 /* rescale probability vector */
  expo->vflag   = 2;                 /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = 2;               /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = RN;                /* N parameters */
  expo->betaVec = Rbeta;             /* beta parameters */
  expo->muVec   = Rmu;               /* mu parameters */
  expo->psiVec  = Rpsi;              /* psi parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */
  expo->est_norm = *estimateNorm;    /* force estimation of matrix norm */
  expo->full_term = 1;

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

  int N = expo->dim;

  expo_alloc_all(expo);
  expo->init_all(expo);

  init_mupsi_const(expo);
  init_wrk(expo);

  int i = 0;

  /********************** expocl ************************/

  int info = 0;
  int one = 1;

  double* p0 = (double*) malloc(expo->dim*sizeof(double));
  double* bak = (double*) malloc(expo->dim*sizeof(double));
  for (i = 0; i < expo->dim; ++i) p0[i] = 1.0;
  memcpy(bak,p0,expo->dim*sizeof(double));

  int n = 5;
  double times[] = { 10.0, 20.0, 30.0, 40.0, 50.0 };
  int ttypes[] = { 0, 0, 0, 0, 0 };

  print_info();

  expo_cl* ecl = expocl_init(expo,CL_DEVICE_TYPE_GPU,128);
  switch (expo->SImodel) {
    case 0:
    case 1:
    default:
      ecl->calc = &expocl_sis_calc;
      break;
    case 2:
      ecl->calc = &expocl_sir_calc;
      break;
  }

  ecl->calc = &expocl_sis_calc;
  expocl_init_kernel(ecl,expo);
  expocl_create_buffers(ecl);
  expocl_create_vector_buffers(ecl);
  expocl_create_index_buffer(ecl,expo);
  printf("\n");

  expocl_copy_mat(ecl,expo);
  expocl_copy_vectors(ecl,expo);
  expocl_copy_index(ecl,expo);

  clReal* b = (clReal*) malloc(expo->dim*sizeof(clReal));
#ifndef DOUBLE
  double_to_clreal(N,p0,b);
  expocl_copy_b(ecl,b);
#else
  expocl_copy_b(ecl,p0);
#endif

  /**************************************************************************/

  double norm_cpu = dnrm2_(&N,p0,&one);
  double norm_gpu = expocl_two_norm_b(ecl);

  /**************************************************************************/

  mach_timespec_t tstart={0,0}, tend={0,0};
  clock_serv_t cclock;

  // Run on CPU
 
  expo->ki = 0;

  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &tstart);
  mach_port_deallocate(mach_task_self(), cclock);

  info = expoTree(n,times,ttypes,p0,expo);

  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &tend);
  mach_port_deallocate(mach_task_self(), cclock);

  printf("     t(CPU) = %fs\n",
         ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - 
         ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
  printf("     # matrix-vector products: %d\n",info);

  norm_cpu = dnrm2_(&N,p0,&one);

  // Run on GPU

  expo->ki = 0;

  memcpy(p0,bak,expo->dim*sizeof(double));

  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &tstart);
  mach_port_deallocate(mach_task_self(), cclock);

  info = clExpoTree(n,times,ttypes,p0,0.0,expo,ecl);

  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &tend);
  mach_port_deallocate(mach_task_self(), cclock);

  printf("     t(GPU) = %fs\n",
         ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - 
         ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
  printf("     # matrix-vector products: %d\n",info);

  expocl_free_kernel(ecl);
  expocl_free(ecl);

  free(b);
  free(p0);
  free(bak);

  /********************** clean up **********************/

  free_wrk(expo);
  expo_free_all(expo);
  expo_type_free(expo);

  return 0;
}

