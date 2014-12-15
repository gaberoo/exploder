#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
using namespace std;

#include <gsl/gsl_math.h>

#include "expoTree.h"
#include "expoTreeSIR.h"

#include "Forest.h"
#include "ExpoTreePars.h"

extern "C" {
  typedef double (*subplx_fun)(int* n, double* x);

  void subplx_(subplx_fun f, int* n, double* tol,
               int* maxnfe, int* mode, double* scale,
               double* x, double* fx, int* nfe,
               double* work, int* iwork, int* iflag);
}

// ===========================================================================

vector<double> K;        /* total population size */
vector<double> beta;     /* infection rate */
vector<double> mu;       /* recovery rate */
vector<double> psi;      /* sampling rate */
vector<double> rShift;   /* time of rate shift */

double rho = 0.0;        /* initial sampling probability */
int root = 0;            /* add a root to the tree */
int n = 0;               /* number of times */
int vflag = 0;           /* vebose flag */
int rescale = 1;
int est_norm = 1;
int extant = 0;
int maxExtant = 0;
int maxEvals = 1000;
double fixedRatio = -1.0;
char jointType = 'x';
char logFun = 'D';
int max_shifts;
int survival = 0;
int nroot = 0;
int model_type = 1;

unsigned char lockVar = 0;

enum Locks {
  LockN    = 0x10,
  LockBeta = 0x08,
  LockMu   = 0x04,
  LockPsi  = 0x02,
  LockRho  = 0x01
};

Forest* T = NULL;

// ===========================================================================

double optim_fun(int* n, double* x) {
  double fx = -INFINITY;

  int i = 0;
  if (! (lockVar & LockN))    K[0]    = x[i++];
  if (! (lockVar & LockBeta)) beta[0] = x[i++];
  if (! (lockVar & LockMu))   mu[0]   = x[i++];
  if (! (lockVar & LockPsi))  psi[0]  = x[i++];
  if (! (lockVar & LockRho))  rho     = x[i++];

  cout << setw(10) << K[0] << " " 
       << setw(10) << beta[0] << " " 
       << setw(10) << mu[0] << " " 
       << setw(10) << psi[0] << " "
       << setw(10) << rho << " :: " << flush;

  if (K[0] < T->maxExtant || beta[0] <= 0.0 || psi[0] < 0.0 || rho <= 0.0 || rho > 1.0) {
    fx = -INFINITY;
  } else {
    if (fixedRatio > 0.0) mu[0] = psi[0]*(1.0/fixedRatio-1.0);
    fx = T->jointLikelihood(K,beta,mu,psi,rho,est_norm,vflag,
                            rescale,survival,model_type);
  }

  cout << " loglik = " << fx << endl;
  if (! gsl_finite(fx)) return INFINITY;
  return -fx;
}

// ===========================================================================

int main(int argc, char* argv[]) {
  int c;
  opterr = 0;
  while ((c = getopt(argc,argv,"n:N:l:vb:u:s:o:r:R:J:S:LE:t:")) != -1) {
    switch (c) {
      case 'n': nroot = atoi(optarg); break;
      case 'N': K.push_back(atof(optarg)); break;
      case 'v': ++vflag; break;
      case 'b': beta.push_back(atof(optarg)); break;
      case 'u': mu.push_back(atof(optarg)); break;
      case 's': psi.push_back(atof(optarg)); break;
      case 'o': rho = atof(optarg); break;
      case 'r': rescale = atoi(optarg); break;
      case 'R': fixedRatio = atof(optarg); break;
      case 'J': jointType = optarg[0]; break;
      case 'S': rShift.push_back(atof(optarg)); break;
      case 'L': logFun = 'L'; break;
      case 'E': est_norm = atoi(optarg); break;
      case 't': lockVar = atoi(optarg); break;
    }
  }

  // =========================================================================

  if (vflag) fprintf(stderr,"Reading times...\n");

  if (optind == argc) {
    fprintf(stderr,"Please provide times file.\n");
    return 0;
  }

  T = new Forest(argc-optind,argv+optind);
  for (size_t i(0); i < rShift.size(); ++i) T->addRateShift(rShift[i]);

  if (vflag) {
    printf("Extant species at time zero = %d\n",T->at(0)->getExtant());
    printf("Max extant species = %d\n",T->at(0)->getMaxExtant());
  }

  // get number of parameters
  int n = 0;
  n += ! (lockVar & LockN);
  n += ! (lockVar & LockBeta);
  n += ! (lockVar & LockMu);
  n += ! (lockVar & LockPsi);
  n += ! (lockVar & LockRho);
  if (vflag) printf("%d free parameters.\n",n);

  double tol = 1e-10;
  int maxnfe = 10000;
  int mode = 0;

  vector<double> x(n);
  int i = 0;
  if (! (lockVar & LockN))    x[i++] = K[0];
  if (! (lockVar & LockBeta)) x[i++] = beta[0];
  if (! (lockVar & LockMu))   x[i++] = mu[0];
  if (! (lockVar & LockPsi))  x[i++] = psi[0];
  if (! (lockVar & LockRho))  x[i++] = rho;

  vector<double> scale(x);

  double fx = -INFINITY;
  int nfe = 0;
  int iflag = 0;

  int nsmax = min(5,n);
  int nsmin = min(2,n);
  int lwork = 2*n + nsmax*(nsmax+4) + 1;
  int liwork = n + int(n/nsmin);

  double* work = (double*) malloc(lwork*sizeof(double));
  int* iwork = (int*) malloc(liwork*sizeof(int));

  subplx_(&optim_fun,&n,&tol,&maxnfe,&mode,scale.data(),x.data(),&fx,&nfe,
          work,iwork,&iflag);

  cout << "iflag = " << iflag << endl;
  cout << "fx = " << fx << endl;
  cout << setw(10) << K[0] << " " 
       << setw(10) << beta[0] << " " 
       << setw(10) << mu[0] << " " 
       << setw(10) << psi[0] << " "
       << setw(10) << rho << endl;

  free(work);
  free(iwork);

  return 0;
}

// ===========================================================================

