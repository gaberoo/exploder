#ifndef __FOREST_H__
#define __FOREST_H__

#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "Tree.h"
#include "expoTree.h"

class Forest : public vector<Tree*> {
  public:
    Forest() : maxExtant(0) {}
    Forest(int n, char** files) : maxExtant(0) { addTrees(n,files); }
    virtual ~Forest() { clear(); }

  // ========================================================================

    inline void addTree(char* fn) { 
      push_back(new Tree(fn)); 
      if (back()->maxExtant > maxExtant) maxExtant = back()->maxExtant;
    }

    inline void addTrees(int n, char** files) {
      for (int i(0); i < n; ++i) addTree(files[i]);
    }

  // ========================================================================

    inline void clear() {
      vector<Tree*>::iterator it(begin());
      for (; it != end(); ++it) delete *it;
      vector<Tree*>::clear();
    }

  // ========================================================================

    inline double likelihood(int treeNum, const vector<double>& vars,
                             int est_norm = 1, int vflag = 0, 
                             int rescale = 1, int survival = 1, 
                             int model = 1) 
    {
      if (treeNum < 0 || treeNum >= (int) size()) return -INFINITY;
      return at(treeNum)->likelihood(vars,est_norm,vflag,rescale,survival,
                                     model);
    }

  // ========================================================================

    inline double jointLikelihood(vector<double>& K, vector<double>& beta, 
                                  vector<double>& mu, vector<double>& psi,
                                  double rho, int est_norm = 1, int vflag = 0,
                                  int rescale = 1, int survival = 1,
                                  int model = 1) 
    {
      double fx(0.0);
      for (size_t i(0); i < size(); ++i) {
        fx += at(i)->likelihood(K,beta,mu,psi,rho,est_norm,
                                vflag,rescale,survival,model);
      }
      return fx;
    }

  inline double jointLikelihood(double K, double beta, double mu, double psi,
      double rho, int est_norm = 1, int vflag = 0, int rescale = 1, int survival = 1) {
    double fx(0.0);
    for (size_t i(0); i < size(); ++i) {
      fx += at(i)->likelihood(K,beta,mu,psi,rho,est_norm,vflag,rescale,survival);
    }
    return fx;
  }

  inline double jointLikelihood(const vector<double>& vars,
      int est_norm = 1, int vflag = 0, int rescale = 1, int survival = 1) {
    double fx(0.0);
    for (size_t i(0); i < size(); ++i) {
      fx += at(i)->likelihood(vars,est_norm,vflag,rescale,survival);
    }
    return fx;
  }

  // ========================================================================

  inline double sumLikelihood(const vector<double>& vars,
      int est_norm = 1, int vflag = 0, int rescale = 1, int survival = 1) {
    double f(0.0);
    double mf(0.0);
    int cnt(0);
    vector<double> fx(size(),0.0);
    for (size_t i(0); i < size(); ++i) {
      fx[i] = at(i)->likelihood(vars,est_norm,vflag,rescale,survival);
      if (isfinite(fx[i])) {
        f += fx[i];
        ++cnt;
      }
    }
    f /= cnt;
    for (size_t i(0); i < size(); ++i) mf += exp(fx[i]-f);
    return log(mf)+f;
  }

  // ========================================================================

  inline double meanLikelihood(const vector<double>& vars,
      int est_norm = 1, int vflag = 0, int rescale = 1, int survival = 1) {
    return sumLikelihood(vars,est_norm,vflag,rescale,survival)-log(size());
  }

  // ========================================================================

  inline void addRateShift(double t) {
    for (size_t i(0); i < size(); ++i) at(i)->addRateShift(t);
  }

  // ========================================================================

  int maxExtant;
};

#endif
