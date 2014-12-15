#include <stdlib.h>

#include "../shared/expo.h"
#include "../shared/expmv.h"

void dgchbv_(int* m, double* t, double* H, int* ldh, double* y, 
             double* wsp, int* iwsp, int* iflag);

int main(int argc, char** argv) {
  double RN[]    = { 100, 200 };
  double Rbeta[] = { 1.0, 2.0 };
  double Rmu[]   = { 0.2, 0.1 };
  double Rpsi[]  = { 0.1, 0.3 };

  double Rki[] = { 10 };
  int estimateNorm[] = { 0 };

  expo_type* expo = (expo_type*) malloc(sizeof(expo_type));

  // allocate parameter structure and initialize
  init_expo(expo);

  expo->N       = (int) ceil(RN[0]); /* total population size */
  expo->K       = RN[0];             /* carrying capacity */
  expo->ki      = *Rki;              /* current number of lineages */
  expo->beta    = Rbeta[0];          /* branching rate */
  expo->mu      = Rmu[0];            /* extinction rate */
  expo->psi     = Rpsi[0];           /* sampling rate */
  expo->SImodel = 1;                 /* Model: DD (SI=1), INF (SI=0) */
  expo->rescale = 1;                 /* rescale probability vector */
  expo->vflag   = 5;                 /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = 2;               /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = RN;                /* N parameters */
  expo->betaVec = Rbeta;             /* beta parameters */
  expo->muVec   = Rmu;               /* mu parameters */
  expo->psiVec  = Rpsi;              /* psi parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */
  expo->est_norm = *estimateNorm;    /* force estimation of matrix norm */

  /* get maximum carrying capacity */
  expo->N_max = max_pop_size(expo);
  int N = expo->N_max+1;

  /* allocate workspace for matrix */
  expo->mat = (double*) malloc(3*N*sizeof(double));

  /* choose lambda function */
  switch (expo->SImodel) {
    case 0:
      expo->lambda = &lambdaInf;
      break;
    default:
      expo->lambda = &lambdaSI;
      break;
  }

  /* allocate workspace for lambda */
  expo->lambdaVec = (double*) malloc(N*sizeof(double));
  init_lambda(expo);

  init_mat(expo);

  double t0[] = { 0.0 };

  int ncol = 1;                  /* EXPMV parameters */

  /* calculate required memory and allocate */
  int memlen = (2*N+(expo->p_max-1)*expo->m_max);
  int wlen   = 2*expo->p_max + (6*ncol+3)*N + ncol + expo->m_max*(expo->p_max-1);
  int iwlen  = 2*N + 4;

  int wrklen = memlen+wlen;
  int iwrklen = iwlen;

  double* wrk = (double*) malloc((memlen+wlen)*sizeof(double));
  int* iwrk = (int*) malloc(iwlen*sizeof(int));

  /* set present time */
  wrk[0] = *t0;
  double wrk_dt = (argc > 1) ? atof(argv[1]) : 0.2;

  double* p0 = wrk;
  double* pT = wrk + N;
  double* tm = wrk + 2*N;
  double* expowrk = wrk + memlen;

  double* p = (double*) calloc(N,sizeof(double));
  int i;
  for (i = expo->ki; i < N; ++i) p[i] = 1.0/(N-expo->ki);

  memcpy(p0,p,N*sizeof(double));
  int one = 1;
  double nrm = dnrm2_(&N,p0,&one);

  for (i = 0; i < N; ++i) {
    // expo->mat[i] = 1.0;
    // expo->mat[N+i] = 0.0;
    // expo->mat[2*N+i] = 0.0;
  }
  expo->shift = expo->trace(expo)/N;  /* shift matrix */

  int info = 0;
  expmv(wrk_dt,N,expo->matvec,expo->norm,expo->trace,p0+expo->offset,
        1,expo->m_max,expo->p_max,tm,1,'d',expo->shift,0,0,expo->vflag,&info,
        expo->est_norm,wrklen-memlen,expowrk,iwrklen,iwrk,expo);

  double* full_mat = (double*) calloc(N*N,sizeof(double));
  double* p2 = (double*) calloc(N,sizeof(double));
  memcpy(p2,p,N*sizeof(double));

  double* wsp = (double*) malloc(2*N*(N+2)*sizeof(double));
  int* iwsp = (int*) malloc(N*sizeof(double));

  for (i = 0; i < N; ++i) {
    full_mat[i*N+i] = expo->mat[N+i] + expo->shift;
    if (i < N-1) full_mat[(i+1)*N+i] = expo->mat[2*N+i];
    if (i > 0) full_mat[(i-1)*N+i] = expo->mat[i];
  }
  dgchbv_(&N,&wrk_dt,full_mat,&N,p2,wsp,iwsp,&info);

  printf("# INFO = %d\n",info);
  for (i = 0; i < N; ++i) {
    printf("%d %g %g %g | %d %d\n",i,p[i],p0[i],p2[i],iwrk[0],iwrk[1]);
  }
 
  free(full_mat);
  free(p2);
  free(wsp);
  free(iwsp);

  /* clean up */
  free(wrk);
  free(iwrk);

  free(expo->mat);
  expo->mat = NULL;

  free(expo->lambdaVec);
  expo->lambdaVec = NULL;

  free(expo);
  expo = NULL;

  return 0;
}

