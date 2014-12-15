#include <stdlib.h>

#include "../shared/expo.h"
#include "../shared/expmv.h"

int main(int argc, char** argv) {
  double RN[]    = { 100,  50 };
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
  expo->N_max = max_pop_size(expo);
  int N = expo->N_max+1;

  double* p  = (double*) calloc(N,sizeof(double));
  double* p0 = (double*) calloc(N,sizeof(double));
  double* pT = (double*) calloc(N,sizeof(double));

  int i;
  for (i = expo->ki; i < N; ++i) p0[i] = 1.0/(N-expo->ki);
  memcpy(p,p0,N*sizeof(double));

  /* allocate workspace for lambda */
  expo->mat = (double*) malloc(3*N*sizeof(double));
  expo->lambdaVec = (double*) malloc(N*sizeof(double));

  init_lambda(expo);
  init_mat(expo);

  (*expo->ft)(p0,pT,expo);
  (*expo->ft)(pT,p0,expo);
  (*expo->ft)(p0,pT,expo);
  for (i = 0; i < N; ++i) {
    printf("%4d %12g %12g %12g\n",
           i,p[i],pT[i],2.0*expo->lambda(i,expo)*((i < N-1) ? p0[i+1] : 0.0));
  }
  
  printf("\n");
  memcpy(p,pT,N*sizeof(double));

  (*expo->fs)(pT,p0,expo);
  (*expo->fs)(p0,pT,expo);
  (*expo->fs)(pT,p0,expo);
  for (i = 0; i < N; ++i) {
    printf("%4d %12g %12g %12g\n",
           i,p[i],p0[i],((i>0) ? expo->psi*pT[i-1] : 0.0));
  }

  free(expo->mat);
  free(expo->lambdaVec);

  free(p);
  free(p0);
  free(pT);

  return 0;
}

