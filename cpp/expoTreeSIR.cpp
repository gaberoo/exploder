#include "expoTreeSIR.h"

bool checkParams(int parVecLen, double* K, double* beta, double* mu,
                 double* psi, double rho, int nroot) 
{
  for (int i = 0; i < parVecLen; ++i) {
    if (K[i] <= 0.0 || beta[i] <= 0.0 || mu[i] < 0.0 || psi[i] < 0.0) {
      return false;
    }
  }

  double maxK = 0.0;
  for (int i = 0; i < parVecLen; ++i) if (maxK < K[i]) maxK = K[i];
  int maxN = (int) ceil(maxK);

  if (nroot < 0 || nroot > maxN) return false;
  if (rho < 0.0 || rho > 1.0) return false;

  return true;
}

//============================================================================

double expoTreeSIR(int parVecLen, double* K, 
                    double* beta, double* mu, double* psi, 
                    double rho, vector<double>& times, vector<int>& ttypes, 
                    int extant, int est_norm, int vflag, int rescale,
                    int nroot)
{
  double fx = -INFINITY;
  int info = 1;

  if (! checkParams(parVecLen,K,beta,mu,psi,rho,nroot)) {
    if (vflag > 0) {
      fprintf(stderr,"Ilegal parameters. Returning inf.\n");
      fprintf(stderr,"N = %g, beta = %g, mu = %g, psi = %g, rho = %g\n",
              K[0],beta[0],mu[0],psi[0],rho);
    }
    return -INFINITY;
  }

  /* count number of sampled lineages */
  int num_sampled = 0;
  for (size_t i(0); i < times.size(); ++i) {
    switch (ttypes[i]) {
      case 0: 
        ++num_sampled; 
        break;
      default: 
        break;
    }
  }

  int ki = 0;
  int nt = times.size();

  double* ptimes = times.data();
  int*    pttypes = ttypes.data();

  double* maxK = max_element(K,K+parVecLen);
  int maxN = ceil(*maxK);

  int dim = sir_index_N(0,maxN,maxN)+1;
  vector<double> p0(dim,0.0);

  double t0 = 0.0;
  double scale = 0.0;

  double sigma = psi[0]/(mu[0]+psi[0]);         /* sampling probability */
  double pR;

  /* set initial value of p */
  if (extant == 0) {
    // p0[0] = 0.0;
    scale = extant*log(sigma);
    int m = 0;
    for (int R = 0; R <= maxN; ++R) {
      pR = (R >= num_sampled) ? pow(1.0-sigma,R-num_sampled) : 0.0;
      p0[m++] = pR;
      for (int I = 1; I <= maxN-R; ++I) {
        p0[m++] = 0.0;
      }
    }
    ki = 0;
    t0 = 0.0;
    ptimes = times.data();
    pttypes = ttypes.data();
    nt = times.size();
  } else {
    cerr << "Not yet implemented for extant lineages !!" << endl;
    return -INFINITY;
  }

#ifdef OPENCL
  info = 2;
  clrExpoTree(K,&ki,beta,mu,psi,&nt,&parVecLen,
              ptimes,pttypes,p0.data(),
              &t0,&info,&est_norm,&vflag,&rescale);
#else
  sir_expotree(K,&ki,beta,mu,psi,&nt,&parVecLen,
               ptimes,pttypes,p0.data(),
               &t0,&info,&est_norm,&vflag,&rescale);
#endif

  if (info > 0) {
    int m = sir_index_N(1,0,maxN);
    fx = p0[m] + scale;
    if (vflag > 0) fprintf(stderr,"ln(p(1,t)) = %20.12e\n",fx);
  } else {
    if (vflag > 0) fprintf(stderr,"rExpoTree returned %d!\n",info);
    return -INFINITY;
  }

  return fx;
}

// ===========================================================================

double expoTreeSIRSurvival(int parVecLen, double* K, double* beta, double* mu, 
                           double* psi, double rho, vector<double>& times, 
                           vector<int>& ttypes, int extant, int est_norm, 
                           int vflag, int rescale, int nroot)
{
  double fx = -INFINITY;

  if (! checkParams(parVecLen,K,beta,mu,psi,rho,nroot)) {
    if (vflag > 0) {
      fprintf(stderr,"Ilegal parameters. Returning inf.\n");
      fprintf(stderr,"N = %g, beta = %g, mu = %g, psi = %g, rho = %g\n",
              K[0],beta[0],mu[0],psi[0],rho);
    }
    return -INFINITY;
  }

  double* maxK = max_element(K,K+parVecLen);
  int maxN = ceil(*maxK);

  int dim = sir_index_N(0,maxN,maxN)+1;
  vector<double> p0(dim,0.0);

  double t0 = 0.0;
  double scale = 0.0;

  int ki = 0;
  int nt = times.size();
  int Ncur = (int) ceil(K[0]);

  p0[0] = 1.0; /* probability that the tree went extinct before the present */
  for (int R(0); R <= Ncur; ++R) {
    for (int I(0); I <= Ncur-R; ++I) {
      int i = sir_index_N(I,R,maxN);
      p0[i] = (rho > 0.0) ? pow(1.-rho,I) : 1.0;
    }
  }

  int info = 0;
  sir_expotree(K,&ki,beta,mu,psi,&nt,&parVecLen,
               times.data(),ttypes.data(),p0.data(),
               &t0,&info,&est_norm,&vflag,&rescale);

  if (info > 0) {
    int m = sir_index_N(1,0,maxN);
    fx = 1.0 - exp(p0[m] + scale);
    if (vflag > 0) fprintf(stderr,"ln(p(1,t)) = %20.12e\n",fx);
  } else {
    if (vflag > 0) fprintf(stderr,"rExpoTree returned %d!\n",info);
    return -INFINITY;
  }

  return (fx > 0) ? log(fx) : -INFINITY;
}


