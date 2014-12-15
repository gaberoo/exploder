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

#include "expoTree.h"

double expoTreeEval(int parVecLen, double* K, 
                    double* beta, double* mu, double* psi, 
                    double rho, vector<double>& times, vector<int>& ttypes, 
                    int extant, int est_norm, int vflag, int rescale,
                    int nroot)
{
  double fx = -INFINITY;
  int info = 1;

  bool goodParams = check_params(parVecLen,K,beta,mu,psi,rho);

  double maxK = 0.0;
  for (int i = 0; i < parVecLen; ++i) if (maxK < K[i]) maxK = K[i];
  int maxN = (int) ceil(maxK);
  if (vflag) printf("Maximum dimension: N = %d\n",maxN);
  if (nroot < 0 || nroot > maxN) goodParams = false;

  if (! goodParams) {
    if (vflag > 0) {
      fprintf(stderr,"Ilegal parameters. Returning inf.\n");
      fprintf(stderr,"N = %g, beta = %g, mu = %g, psi = %g, rho = %g\n",
              K[0],beta[0],mu[0],psi[0],rho);
      return -INFINITY;
    }
  } 
  
  {
    int ki = 0;
    int nt = times.size();

    double* ptimes = times.data();
    int*    pttypes = ttypes.data();

    double* p0 = (double*) malloc((maxN+1)*sizeof(double));
    // vector<double> p0(maxN+1);

    double t0 = 0.0;
    double scale = 0.0;

    // set initial value of p
    if (extant == 0) {
      p0[0] = 0.0;
      if (ttypes[0] == 0) {
        ki = 1;
        for (int m = 1; m <= maxN; ++m) p0[m] = psi[0];
      } else {
        ki = 0;
        for (int m = 0; m <= maxN; ++m) p0[m] = 1.0;
      }
      t0 = times[0];
      ptimes = times.data()+1;
      pttypes = ttypes.data()+1;
      nt = times.size()-1;
    } else {
      ki = extant;
      p0[0] = 0.0;
      scale = extant*log(rho);
      for (int m = 1; m <= maxN; ++m) {
        if (m < extant) p0[m] = 0.0;
        else p0[m] = pow(1.-rho,m-extant);
      }
      ptimes = times.data();
      pttypes = ttypes.data();
      nt = times.size();
    }
#ifdef OPENCL
    clrExpoTree(K,&ki,beta,mu,psi,&nt,&parVecLen,
                ptimes,pttypes,p0,
                &t0,&info,&est_norm,&vflag,&rescale);
#else
    rExpoTree(K,&ki,beta,mu,psi,&nt,&parVecLen,
              ptimes,pttypes,p0,
              &t0,&info,&est_norm,&vflag,&rescale);
#endif
    if (info > 0) {
      fx = p0[nroot+1] + scale;
      if (vflag > 0) fprintf(stderr,"ln(p(1,t)) = %20.12e\n",fx);
    } else {
      if (vflag > 0) fprintf(stderr,"rExpoTree returned %d!\n",info);
      return -INFINITY;
    }

    free(p0);
  }

  return fx;
}

// ===========================================================================

double expoTreeSurvival(vector<double>& K, vector<double>& beta, 
    vector<double>& mu, vector<double>& psi, double rho,
    vector<double>& times, vector<int>& ttypes, 
    int extant, int SI, int vf, int rs) 
{
  double fx = -INFINITY;
  int est_norm = 1;

  bool goodParams(true);
  for (size_t i(0); i < beta.size() && goodParams; ++i) {
    if (beta[i] <= 0.0 || mu[i] < 0.0 || psi[i] < 0.0) goodParams = false;
  }
  if (! goodParams || rho < 0.0 || rho > 1.0) {
    if (vf > 0) fprintf(stderr,"Ilegal parameters. Returning inf.\n");
  } else {
    vector<double>::const_iterator maxK(max_element(K.begin(),K.end()));
    int maxN = (int) ceil(*maxK);
    int ki = 0;
    int nt = times.size();
    int parVecLen = beta.size();
    double t0 = 0.0;
    double* p0 = (double*) malloc((maxN+1)*sizeof(double));
    p0[0] = 1.0; /* probability that the tree went extinct before the present */
    for (int i(1); i <= maxN; ++i) {
      /* probability zero lineages were sampled given i exist */
      p0[i] = (rho > 0.0) ? pow(1.-rho,i) : 1.0;
    }
    nt = times.size();
    rExpoTree(K.data(),&ki,beta.data(),mu.data(),psi.data(),&nt,&parVecLen,
              times.data(),ttypes.data(),p0,&t0,&SI,&est_norm,&vf,&rs);
    fx = p0[1];
    if (vf > 0) fprintf(stderr,"p(no samp) = %20.12e\n",exp(fx));
    fx = 1.-exp(fx);
    free(p0);
    p0 = NULL;
  }

  return (fx > 0) ? log(fx) : -INFINITY;
}

// ===========================================================================

void readTimes(string fn, vector<double>& times, vector<int>& ttypes,
    int& extant, int& maxExtant) {
  ifstream in(fn.c_str());
  if (! in.is_open()) {
    cerr << "Problem opening file '" << fn << "'. Aborting." << endl;
    return;
  }
  string input;
  double x1;
  int    x2;
  while (! getline(in,input).eof()) {
    if (input.length() == 0) continue;
    if (input[0] == '#') continue;
    istringstream istr(input);
    istr >> x1 >> x2;
    times.push_back(x1);
    ttypes.push_back(x2);
  }
  in.close();
  // make sure times are sorted
  vector<size_t> p(times.size());
  gsl_sort_index(p.data(),times.data(),1,times.size());
  // apply sort
  gsl_permute(p.data(),times.data(),1,times.size());
  gsl_permute_int(p.data(),ttypes.data(),1,ttypes.size());
  extant = 0;
  maxExtant = 0;
  int n(times.size());
  for (int i(n-1); i >= 0; --i) {
    switch (ttypes[i]) {
      case 1:
      case 11:
        ++extant;
        break;
      case 0:
      case 2:
      case 5:
      case 10:
        --extant;
        break;
      default:
        break;
    }
    if (maxExtant < extant) maxExtant = extant;
  }
  if (extant < 0) {
    fprintf(stderr,"Invalid tree: More samples than infections.\n");
    abort();
  } 
}

// ===========================================================================

double infTreeSurvival(double beta, double mu, double psi, double rho, double torig) 
{
  if (beta <= 0.0 || mu < 0.0 || psi < 0.0) return -INFINITY;
  double c1;
  double c2;
  double p0;
  c1 = beta-mu-psi;
  c1 = sqrt(c1*c1 + 4*beta*psi);
  c2 = -(beta-mu-psi-2*beta*rho)/c1;
  p0 = exp(-c1*torig)*(1.-c2);
  p0 = beta+mu+psi+c1*(p0-(1.+c2))/(p0+1.+c2);
  p0 = p0/(2.*beta);
  p0 = 1-p0;
  return log((p0 >= 0) ? ((p0 <= 1) ? p0 : 1) : 0);
}

// ===========================================================================

double infTreeEval(double beta, double mu, double psi, double rho,
    const vector<double>& times, const vector<int>& ttypes, 
    int extant, int SI, int survival, int vf) 
{
  if (beta <= 0.0 || mu < 0.0 || psi < 0.0) return -INFINITY;

  double c1;
  double c2;
  double lik;
  double p0;

  c1  = beta-mu-psi;
  c1  = sqrt(c1*c1 + 4*beta*psi);
  c2  = -(beta-mu-psi-2*beta*rho)/c1;
  lik = -log(2*beta);
  if (survival) {
    p0 = exp(-c1*times.back())*(1.-c2);
    p0 = beta+mu+psi+c1*(p0-(1.+c2))/(p0+1.+c2);
    p0 = p0/(2.*beta);
    lik -= log(1.-p0);
  }

  if (extant > 0) lik += extant*log(rho);
  for (size_t i(0); i < times.size(); ++i) {
    if (ttypes[i]) lik += log(2*beta/bdss_q(times[i],c1,c2));
    else lik += log(psi*bdss_q(times[i],c1,c2));
  }

  return lik;
}

// ===========================================================================

int check_params(int n, double* N, double* beta, double* mu, 
                 double* psi, double rho) 
{
  int good_params = 1;
  for (int i = 0; i < n && good_params; ++i) {
    if (N[i]    <= 0.0) good_params = 0;
    if (beta[i] <= 0.0) good_params = 0;
    if (mu[i]   <  0.0) good_params = 0;
    if (psi[i]  <  0.0) good_params = 0;
  }
  if (rho < 0.0 || rho > 1.0) good_params = 0;
  return good_params;
}

