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

#include "cl_expomv.h"

int clExpoTree(int n, double* times, int* ttypes, double* p, 
               double t, expo_type* expo, expo_cl* ecl) 
{
  /***************** SI PARAMETERS *****************/

  double dt    = 0.0;       /* time step */
  double nrm   = 0.0;       /* vector norm */
  double scale;             /* log scaling factor */
  int    info  = 0;         /* info variable */

  /******** INTEGRATE DIFFERENTIAL EQUATION ********/

  int i;

  /* setup vector for device use */
  clReal* b = NULL;

#ifndef DOUBLE
  /* if device is single precision, copy vector */
  b = (clReal*) malloc(expo->dim*sizeof(clReal));
  double_to_clreal(expo->dim,p,b);
  /* copy probability vector to device */
  expocl_copy_b(ecl,b);
#else
  // expocl_copy_b(ecl,p);
  // memcpy(b,p,expo->dim*sizeof(double));
  b = p;
#endif

  expocl_copy_b(ecl,b);            /* copy probability vector to device */
  scale = 0.0;                     /* initialize variables */
  nrm = expocl_two_norm_b(ecl);    /* calculate 2-norm of vector */

  /* rescale probability vector if requested */
  if (expo->rescale && nrm != 1.0) {
    scale += log(nrm);
    nrm = 1./nrm;
    expocl_smult_b(ecl,nrm); // dscal_(&N,&nrm,p0,&one);
    nrm = 1.0;
  }

  /* loop over tree events */
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
      if (expo->vflag > 1) {
        Rprintf("%6d/%02d %8.4f/%8.4f %5d ",i,ttypes[i],dt,times[i],expo->ki);
      }

      /* integrate probability vector over interval */
      info = cl_integrate_interval(dt,expo,ecl,b,&scale,&nrm);

      if (expo->vflag > 1) {
        Rprintf("%10.4e % 11.4e | %5d %5d %5d\n",nrm,scale,info,expo->iwrk[0],expo->iwrk[1]);
      }

      if (info < 0) {
        reset_p(p,expo);
        return info;
      }
    } else {
      if (expo->rescale) {
        nrm = expocl_two_norm_b(ecl); // nrm = dnrm2_(&N,p0,&one);
      }
      if (expo->vflag > 1) {
        Rprintf("%6d/%02d %8.4f/%8.4f %5d %10g % 11.4e\n",
                i,ttypes[i],times[i],dt,expo->ki,nrm,scale);
      }
    }

    /* update initial condition with event information */
    if (i < n-1) {
      cl_update_vector(ttypes[i],expo,ecl);
    }

    nrm = expocl_two_norm_b(ecl); // nrm = dnrm2_(&N,p0,&one);

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
      expocl_smult_b(ecl,nrm); // dscal_(&N,&nrm,p0,&one);
    }

    t = times[i];
  }

#ifndef DOUBLE
  expocl_read_b(ecl,b);
  for (i = 0; i < expo->dim; ++i) p[i] = b[i];
  free(b);
#else
  expocl_read_b(ecl,p);
#endif

  if (expo->vflag > 1) {
    switch (expo->SImodel) {
      case 1:
        Rprintf("p0(%d) = %g, scale = %g\n",1,p[sir_index(1,0,expo)],scale);
        break;

      default:
        Rprintf("p0(%d) = %g, scale = %g\n",1,p[1],scale);
        break;
    }
  }

  log_p(p,expo,scale);

  return info;
}

