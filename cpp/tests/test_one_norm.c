#include <stdlib.h>

#include "../shared/expo.h"
#include <gsl/gsl_statistics.h>

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

  expo->dim = N;

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

  expo->init_all(expo);

  double* colSums = (double*) calloc(expo->dim,sizeof(double));

  int failed = 0;
  int i;
  for (i = 0; i < expo->dim; ++i) {
    double a = (-i*(expo->lambda(i,expo)+expo->psi+expo->mu)-expo->shift);
    double b = (i < expo->dim-1) ? ((i+1)-expo->ki)*expo->mu : 0.0;
    double c = (i > expo->ki) ? ((i-1)+expo->ki)*expo->lambda(i-1,expo) : 0.0;

    double aa = expo->mat[i+expo->dim];
    double bb = (i < expo->dim-1) ? expo->mat[i+1] : 0.0;
    double cc = (i > 0) ? expo->mat[2*expo->dim+i-1] : 0.0;

    double c1 = fabs(a)  + fabs(b)  + fabs(c);
    double c2 = fabs(aa) + fabs(bb) + fabs(cc);

    if (i >= expo->ki && c1 != c2) ++failed;
    colSums[i] = c2;

    printf("%4d | %8.2f %8.2f %8.2f | %8.2f || %8.2f %8.2f %8.2f | %8.2f\n",
        i,aa,bb,cc,c2,a,b,c,c1);
  }

  double one_norm = matFuncOneNormAnalytical(expo);
  double max_val  = gsl_stats_max(colSums,1,expo->dim);

  printf("1-norm = %8.2f\n",one_norm);
  printf("max    = %8.2f\n",max_val);

  if (failed <= 0 && one_norm == max_val) {
    printf("\e[1;32mTest passed!\e[0m\n");
  } else {
    printf("\e[1;31mTest failed: %d failues.\e[0m\n",failed);
  }

  free(colSums);

  /* clean up */
  free(expo->mat);
  expo->mat = NULL;

  free(expo->lambdaVec);
  expo->lambdaVec = NULL;

  free(expo);
  expo = NULL;

  return 0;
}

