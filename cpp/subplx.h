#include <stdlib.h>
#include <stdio.h>

typedef double (*subplx_fun)(int* n, double* x);

void subplx_(subplx_fun f, int* n, double* tol,
             int* maxnfe, int* mode, double* scale,
             double* x, double* fx, int* nfe,
             double* work, int* iwork, int* iflag);
/*
INPUT

  f      - user supplied function f(n,x) to be optimized,
           declared external in calling routine

  n      - problem dimension

  tol    - relative error tolerance for x (tol .ge. 0.)

  maxnfe - maximum number of function evaluations

  mode   - integer mode switch with binary expansion
           (bit 1) (bit 0) :
           bit 0 = 0 : first call to subplx
                 = 1 : continuation of previous call
           bit 1 = 0 : use default options
                 = 1 : user set options

  scale  - scale and initial stepsizes for corresponding
           components of x
           (If scale(1) .lt. 0.,
           abs(scale(1)) is used for all components of x,
           and scale(2),...,scale(n) are not referenced.)

  x      - starting guess for optimum

  work   - double precision work array of dimension .ge.
           2*n + nsmax*(nsmax+4) + 1
           (nsmax is set in subroutine subopt.
           default: nsmax = min(5,n))

  iwork  - integer work array of dimension .ge.
           n + int(n/nsmin)
           (nsmin is set in subroutine subopt.
           default: nsmin = min(2,n))

OUTPUT

  x      - computed optimum

  fx     - value of f at x

  nfe    - number of function evaluations

  iflag  - error flag
           = -2 : invalid input
           = -1 : maxnfe exceeded
           =  0 : tol satisfied
           =  1 : limit of machine precision
           =  2 : fstop reached (fstop usage is determined
                  by values of options minf, nfstop, and
                  irepl. default: f(x) not tested against
                  fstop)
           iflag should not be reset between calls to
           subplx.
*/
