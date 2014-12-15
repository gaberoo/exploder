#include <stdlib.h>

#include "../shared/expo.h"

int main() {

  double RN[]    = { 100, 200 };
  double Rbeta[] = { 1.0, 2.0 };
  double Rmu[]   = { 0.2, 0.1 };
  double Rpsi[]  = { 0.1, 0.3 };

  double Rki[] = { 10 };
  int estimateNorm[] = { 1 };

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

  /* get maximum carrying capacity */
  int j;
  int N = (int) ceil(RN[0]);
  int N2;
  for (j = 0; j < expo->parVecLen; ++j) {
    N2 = (int) ceil(RN[j]);
    if (N < N2) N = N2;
  }
  expo->N_max = N;
  N += 1;

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

  expo->shift = expo->trace(expo)/N;  /* shift matrix */

  init_mat(expo);

  double tmp;
  double tmp2;
  double cs = 0.0;
  double cs2 = 0.0;
  int m = 0;
  for (; m <= expo->N_max; ++m) {
    tmp = matColSum(m,expo);
    if (tmp > cs) cs = tmp;
    tmp2 = fabs(expo->mat[N+m]);
    if (m > 0) tmp2 += fabs(expo->mat[2*N+m-1]);
    if (m < N) tmp2 += fabs(expo->mat[m+1]);
    if (tmp2 > cs2) cs2 = tmp2;
  }

  double maxM = .5*expo->N*(1.+(expo->mu+.5*expo->psi)/expo->beta)
                - .5 - .25*expo->ki;
  int maxMi = rint(maxM);
  if (maxMi < expo->ki) maxMi = expo->ki;
  if (maxMi > expo->N) maxMi = expo->N;

  printf("one-norm (LAPACK)        = %g\n",expo->norm(expo));
  printf("one-norm (max col sum)   = %g\n",cs);
  printf("one-norm (max col sum 2) = %g\n",cs2);
  printf("one-norm (analytical)    = %g\n",matColSum(maxMi,expo));

  printf("inf-norm (LAPACK)        = %g\n",matFuncInfNorm(expo));

  /* clean up */
  free(expo->mat);
  expo->mat = NULL;

  free(expo->lambdaVec);
  expo->lambdaVec = NULL;

  free(expo);
  expo = NULL;

  return 0;
}

