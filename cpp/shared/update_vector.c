#include <stdlib.h>
#include <string.h>
#include "expo.h"
#include "expo_type.h"

int update_vector(int ttype, double* p0, double* pT, expo_type* expo) {
  int info = 0;
  int N = expo->dim;
  int one = 1;
  double nrm = 0.0;
  switch (ttype) {
    case 1:
      // branching event
      (*expo->ft)(p0,pT,expo);
      --expo->ki;
      break;

    case 0:
      // sampling-removal event
      (*expo->fs)(p0,pT,expo);
      ++expo->ki;
      break;

    case 2:
      // sample extinct lineage
      ++expo->ki;
    case 3:
      // sampling extant lineage
      (*expo->fs2)(p0,pT,expo);
      break;

    case 4:
      // add lineages without altering likelihood
      memcpy(pT,p0,expo->dim*sizeof(double));
      pT[expo->ki] = 0.0;
      if (++expo->ki > expo->N) {
        if (expo->vflag > 1) {
          Rprintf("Increasing lineages above current maximum (N = %d). Aborting.\n",expo->N);
        }
        info = -5;
      }
      break;

    case 5:
      // extinction by sampling
      memcpy(pT,p0,N*sizeof(double));
      break;

    case 10:
      /* shift vector upwards without altering likelihood */
      memcpy(pT+1,p0,(expo->dim-1)*sizeof(double));
      pT[0] = 0.0;
      ++expo->ki;
      break;

    case 11:
      /* shift vector downwards without altering likelihood */
      memcpy(pT,p0+1,(expo->dim-1)*sizeof(double));
      pT[expo->dim-1] = 0.0;
      --expo->ki;
      break;

    /**** RATE SHIFTS ****/
    case 20: /* immediate extinction */
    case 21:
    case 22:
      if (++expo->curPar < expo->parVecLen) {
        if (expo->vflag > 1) {
          Rprintf("Rate shift!\n");
          if (! expo->user_funcs) {
            Rprintf("  OLD: N = %f, beta = %f, mu = %f, psi = %f, ki = %d\n",
                    expo->K,expo->beta,expo->mu,expo->psi,expo->ki);
          }
          Rprintf("    +++ %g %g %g %g\n",p0[expo->ki-2],p0[expo->ki-1],p0[expo->ki],p0[expo->ki+1]);
        }

        if (! expo->user_funcs) {
          expo->K    = expo->NVec[expo->curPar];
          expo->N    = (int) ceil(expo->K);
          expo->beta = expo->betaVec[expo->curPar];
          expo->mu   = expo->muVec[expo->curPar];
          expo->psi  = expo->psiVec[expo->curPar];
        } else {
          expo->lambdaVec += expo->N_max+1;
          expo->muFun += expo->N_max+1;
          expo->psiFun += expo->N_max+1;
        }

        if (expo->vflag > 1) {
          nrm = dnrm2_(&N,p0,&one);
          if (! expo->user_funcs) {
            Rprintf("  NEW: N = %f, beta = %f, mu = %f, psi = %f. |v| = %f.\n",
                    expo->K,expo->beta,expo->mu,expo->psi,nrm);
          }
          Rprintf("    +++ %g %g %g %g\n",p0[expo->ki-2],p0[expo->ki-1],p0[expo->ki],p0[expo->ki+1]);
        }
      } else {
        if (expo->vflag > 1) {
          Rprintf("Not enough parameters supplied to shift rates. Ignoring rate shift!\n");
        }
      }

      memcpy(pT,p0,N*sizeof(double));
      break;

    case 99:
    default:
      memcpy(pT,p0,N*sizeof(double));
      break;
  }

  return info;
}

