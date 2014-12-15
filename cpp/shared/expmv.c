/*
 * Copyright (c) 2012-2014, Gabriel Leventhal, ETH Zurich
 * All rights reserved.
 *
 * Original MATLAB implementation:
 * Copyright (c) 2010, Nick Higham and Awad Al-Mohy
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

#include "expmv.h"

inline void exchange(double** x, double** y) {
  double* tmp = *x;
  *x = *y;
  *y = tmp;
}

inline void exchange_cpy(int n, double* x, double* y) {
  double* tmp = (double*) malloc(n*sizeof(double));
  memcpy(tmp,y,n*sizeof(double));
  memcpy(y,x,n*sizeof(double));
  memcpy(x,tmp,n*sizeof(double));
}

int expmv(double t, double* b, int recalcm, expo_type* expo)
{
  int n = expo->dim;
  int ncol = 1;

  int mv = 0;           /* number of matrix-vector products computed */
  int mvd = 0;

  int i = 0;
  int k = 0;

  double tt;
  int tcol;         /* size of Taylor approximation */

  int wlen;
  int nalen;

  int   nncol = n*ncol;
  int    ione = 1;
  double done = 1.0;

  /* check supplied memory */
  wlen = 2*expo->p_max + 2*nncol + expo->m_max*(expo->p_max-1);
  nalen = 3*n + (4*n+1)*ncol;
  if (expo->wlen < wlen+nalen || expo->iwlen < 2*n + 4) {
#ifdef DEBUG
    if (expo->vflag) {
      Rprintf("Not enough workspace supplied!\n");
      Rprintf("  Required 'double': %d\n",wlen+nalen);
      Rprintf("  Supplied         : %d\n",wrklen);
      Rprintf("  Required 'int'   : %d\n",2*n+4);
      Rprintf("  Supplied         : %d\n",iwrklen);
    }
#endif
    return -2;
  }

  /* need 2N offset on workspace because that's where p0 and pT are stored! */
  double* talpha = expo->wrk + 2*expo->dim;          /* length = p_max */
  double* teta   = talpha + expo->p_max;             /* length = p_max */
  double* C      = teta + expo->p_max;               /* length = m_max*(p_max-1) */
  double* b1     = C + expo->m_max*(expo->p_max-1);  /* length = n*ncol */
  double* b2     = b1 + nncol;                       /* length = n*ncol */
  double* nawrk  = b2 + nncol;                       /* length = nalen = (4*n+1)*ncol + 3*n > n*ncol */

  /* assign tolerance */

#ifdef DEBUG
  clock_t tx = clock();
  clock_t t_infnorm = 0;
  clock_t t_matvec  = 0;
  size_t  c_infnorm = 0;
  size_t  c_matvec  = 0;
#endif

  // get required Taylor truncation (if not set)
  if (recalcm) {
#ifdef DEBUG
    if (expo->vflag > 2) Rprintf("Calculating required Taylor truncation...");
    tx = clock();
#endif
    tt = 1.;
    select_taylor_degree(t,b,expo,expo->est_norm,talpha,teta,nawrk,expo->iwrk);
#ifdef DEBUG
    int j;
    tx = clock() - tx;
    fprintf(stderr,"done [t = %fs]\n",(float) tx/CLOCKS_PER_SEC);

    if (expo->vflag > 4) {
      for (i = 0; i < m_max; ++i) {
        for (j = 0; j < p_max-1; ++j) {
          Rprintf("%16.8e ",expo->tm[j*m_max+i]);
        }
        Rprintf("\n");
      }
    }
#endif
    mv = mvd;
  } else {
    tt = t; 
    mv = 0; 
    mvd = 0;
  }

  double s = 1.0;      /* cost per column */

  if (t == 0.0) { 
    tcol = 0; 
  } else {
#ifdef DEBUG
    if (expo->vflag > 3) Rprintf("Calculating cost...");
    tx = clock();
#endif
    taylor_length(tt,expo->tm,expo,&s,&tcol);
#ifdef DEBUG
    tx = clock() - tx;
    if (expo->vflag > 3) Rprintf("done [t = %fs]\n",(float) tx/CLOCKS_PER_SEC);
#endif
  }

  if (tcol == 0) {
    if (expo->vflag) {
      Rprintf("Cannot calculate matrix exponential (under-/overflow?).\n");
      Rprintf("Returned results may be gibberish.\n");
    }
    return -1;
  }

  /* scaling factor eta */
  double eta = (expo->shift != 0.0) ? exp(t*expo->shift/s) : 1.0;

  memcpy(b1,b,nncol*sizeof(double));

#ifdef DEBUG
  if (expo->vflag > 2) {
    Rprintf("m = %2d, s = %g, ||b|| = %g\n", 
            tcol, s, inf_norm(n,ncol,b));
  }
#endif

  double c1 = 0.0;
  double c2 = 0.0;
  double bnorm = 0.0;

  int ss;

  for (ss = 1; ss <= s; ++ss) {
    /* get norm of input vector/matrix */
#ifdef DEBUG
    tx = clock();
    if (expo->vflag > 2) Rprintf("Calculating norm...");
#endif
    c1 = inf_norm(n,ncol,b1);
#ifdef DEBUG
    if (expo->vflag > 2) Rprintf("done. |x| = %8.2f.\n",c1);
    t_infnorm += clock() - tx;
    ++c_infnorm;
#endif

    for (k = 1; k <= tcol; ++k) {
      /* apply matrix to vector, but scale by t/(s*k) */
#ifdef DEBUG
      tx = clock();
#endif
      expo->matvec('n',n,ncol,t/(s*k),b1,b2,expo);

#ifdef DEBUG
      t_matvec += clock() - tx;
      ++c_matvec;
#endif
      ++mv;  /* increment matrix-vector product counter */

      /* add b2 onto b */
      daxpy_(&nncol,&done,b2,&ione,b,&ione);

      /* check whether convergence has been reached before full truncation */
      if (! expo->full_term) {
#ifdef DEBUG
        tx = clock();
#endif
        /* norm of output vector */
        c2 = inf_norm(n,ncol,b2);
#ifdef DEBUG
        t_infnorm += clock() - tx;
        ++c_infnorm;

        // if (expo->vflag > 2) Rprintf("k=%3d: |b| = %9.2e, c1 = %9.2e, c2 = %9.2e.\n",k,bnorm,c1,c2);
#endif

#ifdef DEBUG
        tx = clock();
#endif
        /* get new infinite norm of the new b */
        bnorm = inf_norm(n,ncol,b);
#ifdef DEBUG
        t_infnorm += clock() - tx;
        ++c_infnorm;
        if (expo->vflag > 3) Rprintf("  k=%3d: |b| = %9.2e, c1 = %9.2e, c2 = %9.2e.\n",k,bnorm,c1,c2);
#endif

        if (c1+c2 <= expo->prec*bnorm) break;

        c1 = c2;
      }

      /* exchange b1 and b2 */
      exchange(&b1,&b2);
    }

#ifdef DEBUG
    if (expo->vflag > 2) {
      Rprintf("m = %3d: %9.2e %9.2e %9.2e %9.2e\n",k,bnorm,c1,c2,(c1+c2)/bnorm);
    }
#endif

    for (i = 0; i < nncol; ++i) b[i] *= eta;

    if (ss < s) memcpy(b1,b,nncol*sizeof(double));
  }

#ifdef DEBUG
  float xt = (float) t_infnorm/CLOCKS_PER_SEC;
  fprintf(stderr,"time in inf-norm = %fs. Averge = %f (%lu).\n",
          xt,xt/c_infnorm,c_infnorm);

  xt = (float) t_matvec/CLOCKS_PER_SEC;
  fprintf(stderr,"time in mat-vec  = %fs. Averge = %f (%lu).\n",
          xt,xt/c_matvec,c_matvec);
#endif

  expo->iwrk[0] = tcol;
  expo->iwrk[1] = k;
  expo->iwrk[2] = s;

  return mv;
}

/****************************************************************************/
 
double find_absmax_array(int i0, int len, const double* x) {
  if (len > 1) {
    int i3 = i0+len/2;
    return fmax(find_absmax_array(i0,i3,x),find_absmax_array(i3+1,i0+len,x));
  } else if (len == 1) {
    return fmax(fabs(x[i0]),fabs(x[i0+len]));
  } else {
    return fabs(x[i0]);
  }
}

/****************************************************************************/

double inf_norm(int n1, int n2, double* A) {
  double c = 0.0;
#if defined(UNROLL)
  /* use unrolled loops */
  int i;
  double rowsum[4];
  int m = n1 % 4;
  if (m != 0) {
    for (i = 0; i < m; ++i) {
      rowsum[i] = fabs(A[i]);
      if (rowsum[i] > c) c = rowsum[i];
    }
  }
  if (n1 < 4) return c;
  int j1, j2;
  for (i = m; i < n1; i += 4) {
    rowsum[0] = fabs(A[i]);
    rowsum[1] = fabs(A[i+1]);
    rowsum[2] = fabs(A[i+2]);
    rowsum[3] = fabs(A[i+3]);
    j1 = (rowsum[0] > rowsum[1]) ? 0 : 1;
    j2 = (rowsum[2] > rowsum[3]) ? 2 : 3;
    if (rowsum[j2] > rowsum[j1]) j1 = j2;
    if (rowsum[j1] > c) c = rowsum[j1];
  }
#elif defined(IDAMAX)
  int one = 1;
  int row = idamax_(&n1,A,&one);
  c = A[row];
#elif defined(DLANGE)
  char infNormChar = 'i';
  c = dlange_(&infNormChar,&n,&ncol,b1,&n,nawrk);
#else
  double rowsum;
  int i;
  for (i = 0; i < n1; ++i) {
    // rowsum = 0.0;
    // int j;
    // for (j = 0; j < n2; ++j) rowsum += fabs(A[j*n1+i]);
    rowsum = fabs(A[i]);
    if (rowsum > c) c = rowsum;
  }
#endif
  return c;
}


