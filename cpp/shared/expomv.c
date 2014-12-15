/*
 * Copyright (c) 2012-2014, Gabriel Leventhal, ETH Zurich
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer
 *     in the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of the ETH Zurich nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "expomv.h"

int expoTree(int n, double* times, int* ttypes, 
              double* p, expo_type* expo) 
{
  /***************** SI PARAMETERS *****************/

  double t = 0.0;           /* total time */
  double dt = 0.0;          /* time step */
  double nrm = 1.0;         /* vector norm */
  double scale = 0.0;       /* log scaling factor */
  int one = 1;              /* integer 1 (FORTRAN corpse) */
  int info = 0;

  /********************* EXPMV *********************/

  int N = expo->dim;        /* maximum dimension */

  /* split up calculation when precision is bad */
//  int    wrk_steps     = 1;
//  double wrk_dt        = 0.0;
//  int    wrk_info      = 0;

  /* initialize variables */
  t = expo->wrk[0];

  /* TODO: too much workspace supplied!! */
  /* set pointers for easier access to variables */
  double* p0 = expo->wrk;
  double* pT = expo->wrk + N;

  /******** INTEGRATE DIFFERENTIAL EQUATION ********/

  /* copy probability vector and calculate norm */
  memcpy(p0,p,N*sizeof(double));
  nrm = dnrm2_(&N,p0,&one);

  /* rescale probability vector if requested */
  if (expo->rescale && nrm != 1.0) {
    scale += log(nrm);
    nrm = 1./nrm;
    dscal_(&N,&nrm,p0,&one);
    nrm = dnrm2_(&N,p0,&one);
  }

  /* loop over tree events */
  int i;
  // int j;

  for (i = 0; i < n; ++i) {
    /* get time interval between events */
    dt = times[i] - t;

    if (dt < 0.0) {
      Rprintf("Negative dt (%f) at time step %d! Aborting.\n",dt,i);
      Rprintf("Numer of time points: n = %d\n",n);
      Rprintf("t(0)   = %8.4e\n",times[0]);
      Rprintf("t(i)   = %8.4e\n",t);
      Rprintf("t(i+1) = %8.4e\n",times[i]);
      Rprintf("dt     = %8.4e\n",dt);
      reset_p(p,expo);
      return -1;
    }

    /* don't do anything for dt = 0 */
    if (dt > 0.0) {
      /* save a backup of the vector */
      memcpy(pT,p0,N*sizeof(double));

      if (expo->vflag > 1) {
        Rprintf("%6d/%02d %8.4f/%8.4f %5d ",i,ttypes[i],dt,times[i],expo->ki);
      }
//      if (expo->vflag > 1) {
//        if (wrk_steps == 1) {
//        } else {
//          Rprintf("%8d %8f %5d %8.4e % 8.4e %2d %2d\n",
//              k,wrk_dt,expo->ki,nrm,scale+wrk_scale,iwrk[0],iwrk[1]);
//        }
//      }

      info = integrate_interval(dt,expo,p0,&scale,&nrm);

      if (expo->vflag > 1) {
        Rprintf("%10.4e % 11.4e | %5d %5d %5d\n",
                nrm,scale,info,expo->iwrk[0],expo->iwrk[1]);
      }

      if (info < 0 || expo->iwrk[0] < 0) {
        reset_p(p,expo);
        expo->iwrk[0] = info;
        return expo->iwrk[0];
      }
    } else {
      nrm = dnrm2_(&N,p0,&one);
      if (expo->vflag > 1) {
        Rprintf("%6d/%02d %8.4f/%8.4f %5d %10g\n",
                i,ttypes[i],times[i],dt,expo->ki,nrm);
      }
    }

    // update initial condition with event information
    if (i < n-1) {
      info = update_vector(ttypes[i],p0,pT,expo);
      if (expo->iwrk[0] < 0) {
        reset_p(p,expo);
        return info;
      }
      memcpy(p0,pT,N*sizeof(double));

      if (expo->p_save != NULL) memcpy(expo->p_save+i*N,pT,N*sizeof(double));
    }

    nrm = dnrm2_(&N,p0,&one);

    if (nrm <= 0.0) {
      if (expo->vflag > 2) {
        Rprintf("Vector norm zero. Illegal event.\n");
      }
      reset_p(p,expo);
      return -6;
    }

    /* rescale likelihoods for numerical reasons */
    if (expo->rescale) {
      if (nrm > 1e20) {
        Rprintf("Problem with 2-norm of vector in rescaling!\n");
        reset_p(p,expo);
        return -7;
      }
      scale += log(nrm);
      nrm = 1./nrm;
      dscal_(&N,&nrm,p0,&one);
    }

    t = times[i];
  }

  if (expo->vflag > 1) {
    Rprintf("p0(%d) = %g, scale = %g\n",1,p0[1],scale);
  }

  memcpy(p,p0,N*sizeof(double));
  log_p(p,expo,scale);

  
  return 1;
}

/************************************************************/ 

void log_p(double* p, const expo_type* expo, double scale) {
  int j;
  for (j = 0; j < expo->N_max+1; ++j) {
    p[j] = log(p[j]) + scale;
  }

}

/************************************************************/ 

void reset_p(double* p, const expo_type* expo) {
  int N = expo->N_max+1;
  int j;
  for (j = 0; j < N; ++j) p[j] = R_NegInf;
}

/************************************************************/ 

void printSaved(int n, const expo_type* expo) {
  size_t i, j;
  size_t N = expo->N_max+1;
  for (i = 0; i < n; ++i) {
    for (j = 0; j < N; ++j) {
      printf("%12g ",expo->p_save[i*N+j]);
    }
    printf("\n");
  }
}

