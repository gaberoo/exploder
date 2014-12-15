#include <stdlib.h>

#include "../shared/expo.h"
#include "../shared/expotree_ext.h"

int main() {
  double RN[]    = { 10 };
  double Rki[] = { 1 };
  int estimateNorm[] = { 1 };

  expo_type* expo = expo_type_alloc();

  expo->N       = (int) ceil(RN[0]); /* total population size */
  expo->K       = RN[0];             /* carrying capacity */
  expo->ki      = *Rki;              /* current number of lineages */
  expo->SImodel = 1;                 /* Model: DD (SI=1), INF (SI=0) */
  expo->rescale = 1;                 /* rescale probability vector */
  expo->vflag   = 0;                 /* verbosity level */
  expo->cutoff  = 0.0;               /* precision for zero */

  expo->parVecLen = 1;               /* number of parameter sets */
  expo->curPar  = 0;                 /* current parameter set */
  expo->NVec    = RN;                /* N parameters */

  expo->offset = 0;                  /* don't calculate for zero entries */
  expo->est_norm = *estimateNorm;    /* force estimation of matrix norm */

  expo->user_funcs = 1;              /* use user-supplied functions */
  expo->trace = &matFuncTrace_ext;   /* matrix trace */
  expo->fs = &sampleEvent_ext;       /* sampling event */
  expo->init_all = &init_mat_ext;    /* only the matrix needs to be init'ed */

  /* get maximum carrying capacity */
  expo->N_max = max_pop_size(expo);
  int N = expo->N_max+1;

  /* allocate workspace for matrix */
  expo->mat = (double*) malloc(3*N*sizeof(double));

  double lambda[] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
  double mu[]     = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
  double psi[]    = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

  /* allocate workspace for lambda */
  expo->lambdaVec = lambda;
  expo->muFun     = mu;
  expo->psiFun    = psi;

  expo->init_all(expo);

  int failed = 0;
  int i;
  for (i = 0; i < N; ++i) {
    double a = (-i*(lambda[i]+psi[i]+mu[i])-expo->shift);
    double b = (i < expo->N) ? (i+expo->ki)*lambda[i] : 0.0;
    double c = (i > expo->ki) ? (i-expo->ki)*mu[i] : 0.0;
    printf("%4d | %8.2f %8.2f %8.2f | %8.2f %8.2f %8.2f\n",
        i,expo->mat[i],expo->mat[N+i],expo->mat[2*N+i],a,b,c);
    if (i >= expo->ki) {
      failed += expo->mat[i] != c;
      failed += expo->mat[N+i] != a;
      failed += expo->mat[2*N+i] != b;
    }
  }

  if (failed <= 0) {
    printf("\e[1;32mTest passed!\e[0m\n");
  } else {
    printf("\e[1;31mTest failed: %d failues.\e[0m\n",failed);
  }

  /* clean up */
  free(expo->mat); expo->mat = NULL;

  expo_type_free(expo);

  return 0;
}

