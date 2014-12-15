#include <stdlib.h>
#include <stdio.h>

#include "../shared/expo.h"
#include "../shared/expmv.h"

#include "expocl.h"
#include "cl_expomv.h"

// #include "tests/start_vec.h"
#include <time.h>
#include <mach/clock.h>
#include <mach/mach.h>

void dgemv_(char* TRANS, int* M, int* N, double* ALPHA, double* A, 
            int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);

inline void inc_clock(clock_t* t, clock_t* last_t) {
  *t += clock() - *last_t;
  *last_t = clock();
}

float snrm2_(int* n, float* x, int* incx);

int main(int argc, char** argv) {
  double RN[]    = { 1000 };
  double Rbeta[] = { 100.0 };
  double Rmu[]   = { 10.0 };
  double Rpsi[]  = { 10.0 };
  double Rki[] = { 1 };
  int estimateNorm[] = { 0 };

  Rbeta[0] /= RN[0];
  Rmu[0]   /= RN[0];
  Rpsi[0]  /= RN[0];

  expo_type* expo = expo_type_alloc();
  init_expo(expo);

  expo->N       = (int) ceil(RN[0]); /* total population size */
  expo->K       = RN[0];             /* carrying capacity */
  expo->ki      = *Rki;              /* current number of lineages */
  expo->beta    = Rbeta[0];          /* branching rate */
  expo->mu      = Rmu[0];            /* extinction rate */
  expo->psi     = Rpsi[0];           /* sampling rate */
  expo->SImodel = (argc > 1) ? atoi(argv[1]) : 1; /* Model: DD (SI=1), INF (SI=0) */
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

  int N = expo->dim;

  expo_alloc_all(expo);
  expo->init_all(expo);
  init_mupsi_const(expo);
  init_wrk(expo);

  int i = 0;
  char trans = 't';

  int ncol = 1;

  /* set present time */
  expo->wrk[0] = 0.0;

  double* p0        = (double*) calloc(expo->dim,sizeof(double));
  double* start_vec = (double*) calloc(expo->dim,sizeof(double));

  for (i = 0; i < expo->dim; ++i) {
    // start_vec[i] = ((double) rand())/((double) RAND_MAX);
    start_vec[i] = 1.0;
  }
  memcpy(p0,start_vec,expo->dim*sizeof(double));

  /******************** full matrix *********************/

  int m;
  double* A = NULL;
  if (expo->dim < 10000) {
    printf("Allocating A matrix (dim=%d)...",expo->dim);
    A = (double*) calloc(expo->dim*expo->dim,sizeof(double));
    if (expo->SImodel == 2) {
      for (m = 0; m < expo->dim; ++m) {
        int I = expo->lookup[m];
        int R = expo->lookup[m+expo->dim];
        int a = sir_index(I-1,R+1,expo);
        int b = sir_index(I+1, R ,expo);
        A[m*expo->dim+m] = expo->mat[m];
        if (a >= 0) A[a*expo->dim+m] = expo->mat[m+expo->dim];
        if (b >= 0) A[b*expo->dim+m] = expo->mat[m+2*expo->dim];
      }
    } else {
      for (m = 0; m < expo->dim; ++m) {
        if (m > 0) A[(m-1)*expo->dim+m] = expo->mat[m];
        A[m*expo->dim+m] = expo->mat[m+expo->dim];
        if (m < expo->dim-1) A[(m+1)*expo->dim+m] = expo->mat[m+2*expo->dim];
      }
    }
    printf("done.\n");
  }

  /********************** expocl ************************/

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

  expocl_init_kernel(ecl,expo);

  expocl_create_buffers(ecl);
  expocl_create_index_buffer(ecl,expo);
  expocl_create_vector_buffers(ecl);

  expocl_copy_mat(ecl,expo);
  expocl_copy_vectors(ecl,expo);
  expocl_copy_index(ecl,expo);

  printf("\n");

  clReal* b = (clReal*) malloc(expo->dim*sizeof(clReal));
  double_to_clreal(N,start_vec,b);

  expocl_copy_x(ecl,b);

  /********************** matvec ************************/

  double cpu_norm = inf_norm(N,ncol,p0);
  double gpu_norm = expocl_inf_norm_x(ecl);
  printf("CPU = %f :: GPU = %f\n",cpu_norm,gpu_norm);

  double* cpu_out = (double*) calloc(expo->dim,sizeof(double));
  double* gpu_out = (double*) calloc(expo->dim,sizeof(double));
  double* blas_out = (double*) calloc(expo->dim,sizeof(double));

  double dt = 10.0;

  mach_timespec_t tstart={0,0}, tend={0,0};
  clock_serv_t cclock;

  /***** CPU *****/

  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &tstart);
  mach_port_deallocate(mach_task_self(), cclock);

  expo->matvec(trans,0,0,dt,p0,cpu_out,expo);

  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &tend);
  mach_port_deallocate(mach_task_self(), cclock);

  printf("     t(CPU) = %fs\n",
         ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - 
         ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
   
  /***** GPU *****/

  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &tstart);
  mach_port_deallocate(mach_task_self(), cclock);

  ecl->calc(ecl,expo,trans,dt);
  expocl_copy_xy(ecl);
  expocl_read_x(ecl,b);

  host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
  clock_get_time(cclock, &tend);
  mach_port_deallocate(mach_task_self(), cclock);

  printf("     t(GPU) = %fs\n",
         ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - 
         ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));

  clreal_to_double(N,b,gpu_out);

  /***** blas ******/

  int one = 1;
  double dzero = 0.0;
  if (A != NULL) {
    dgemv_(&trans,&expo->dim,&expo->dim,&dt,A,&expo->dim,start_vec,&one,&dzero,blas_out,&one);
  }

  for (i = 0; i < expo->dim; ++i) {
    // int I = expo->lookup[i];
    // int R = expo->lookup[i+expo->dim];
    if (cpu_out[i] != gpu_out[i]) {
      printf("%4d %30.26f %30.26f %30.26f\n",i,cpu_out[i],gpu_out[i],blas_out[i]);
    }
  }

  cpu_norm = inf_norm(N,ncol,cpu_out);
  gpu_norm = expocl_inf_norm_x(ecl);
  double gpn2 = inf_norm(N,ncol,gpu_out);
  double nrm_blas = inf_norm(expo->dim,ncol,blas_out);

//  gpu_norm = expocl_two_norm_x(ecl);
//  double nrm_blas = dnrm2_(&expo->dim,blas_out,&one);

  printf("CPU = %g :: GPU = %g (%g) :: BLAS = %g\n",cpu_norm,gpu_norm,gpn2,nrm_blas);

  expocl_free_kernel(ecl);
  expocl_free(ecl);

  free(b);
  free(p0);
  free(start_vec);
  if (A != NULL) free(A);

  /********************** clean up **********************/

  free(expo->mat);       expo->mat = NULL;
  free(expo->mat_i);     expo->mat_i = NULL;
  free(expo->lookup);    expo->lookup = NULL;
  free(expo->lambdaVec); expo->lambdaVec = NULL;
  free(expo->muFun);     expo->muFun = NULL;
  free(expo->psiFun);    expo->psiFun = NULL;

  free_wrk(expo);
  expo_type_free(expo);

  return 0;
}

