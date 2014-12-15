#ifndef __TREE_H__
#define __TREE_H__

#include "expoTree.h"
#include "expoTreeSIR.h"
// #include "TreeNode.h"

#include <vector>
#include <string>
using namespace std;

class Tree {
public:
  /* CONSTRUCTORS-DESTRUCTORS */

  // 1. Default
  Tree() : extant(0), maxExtant(0), nroot(0), fn(""), is_rev(0), numBranches(0) {}

  // 2. Read tree from file
  Tree(const char* fn) : nroot(0), is_rev(0) { read(string(fn)); }
  Tree(string fn) : nroot(0), is_rev(0) { read(fn); }

  // 3. Copy constructor
  Tree(const Tree& T) 
    : times(T.times), ttypes(T.ttypes), extant(T.extant),
      maxExtant(T.maxExtant), nroot(T.nroot), is_rev(T.is_rev),
      numBranches(T.numBranches) //, nodes(T.nodes)
  {
    ids.resize(T.ids.size());
    for (size_t i(0); i < T.ids.size(); ++i) ids[i] = T.ids[i];
  }

  // Destructor
  virtual ~Tree() {}

  // ========================================================================
  // METHODS
  
  // read tree from file
  inline void read(string filename) { 
    fn = filename;
    // readTimes(fn,times,ttypes,extant,maxExtant); 
    readFromFile(fn);
  }
  void readFromFile(string fn);

  // get pointers
  inline double* ti() { return times.data(); }
  inline int* tt() { return ttypes.data(); }

  inline size_t size() const { return times.size(); }

  // set methods
  inline void set_nroot(int n) { nroot = n; }

  // get methods
  inline double time(size_t i) const { return times.at(i); }
  inline int ttype(size_t i) const { return ttypes.at(i); }
  inline int getMaxExtant() const { return maxExtant; }
  inline int getExtant() const { return extant; }

  inline size_t countTypes(int type) const {
    size_t c = 0;
    for (size_t i = 0; i < ttypes.size(); ++i) {
      if (ttype(i) == type) ++c;
    }
    return c;
  }

  inline int num_branches() const { return numBranches; }

  inline int max_id() const {
    int max = 0;
    for (size_t i = 0; i < ids.size(); ++i) {
      if (ids[i].size() > 0) {
        if (ids[i][0] > max) max = ids[i][0];
      }
    }
    return max;
  }

  // ========================================================================
  // OPERATORS

  const Tree& operator=(const Tree& t) {
    times = t.times;
    ttypes = t.ttypes;
    extant = t.extant;
    maxExtant = t.maxExtant;
    fn = t.fn;
    return *this;
  }

  // ========================================================================

  inline double maxTime() const { return times.back(); }
  void reverse();

  // ========================================================================

  // calculate survival probability of the tree
  double survival(vector<double>& K, vector<double>& beta, 
      vector<double>& mu, vector<double>& psi, double rho,
      int est_norm = 1, int vflag = 0, int rescale = 1, int nroot = 0, 
      int model_type = 1);

  inline double survival(double K, double beta, double mu, 
      double psi, double rho, int est_norm = 1, int vflag = 0, 
      int rescale = 1) 
  {
    double fx(-INFINITY);
    if (K < 1.0) 
      fx = infTreeSurvival(beta,mu,psi,rho,times.back());
    else 
      fx = expoTreeSurvival(K,beta,mu,psi,rho,times.back(),
                            extant,est_norm,vflag,rescale);
    return fx;
  }

  // ========================================================================

  double likelihood(vector<double>& K, vector<double>& beta, 
                    vector<double>& mu, vector<double>& psi, 
                    double rho, int est_norm = 1, int vflag = 0, 
                    int rescale = 1, int surv = 1, int model = 1);

  inline double likelihood(double K, double beta, double mu, double psi, 
                           double rho, int est_norm = 1, int vflag = 0, 
                           int rescale = 1, int surv = 1, int model = 1)
  {
    vector<double> _K(1,K);
    vector<double> _b(1,beta);
    vector<double> _m(1,mu);
    vector<double> _p(1,psi);
    return likelihood(_K,_b,_m,_p,rho,est_norm,vflag,rescale,surv,model);
  }

  inline double likelihood(const vector<double>& vars,
                           int est_norm = 0, int vflag = 0, int rescale = 1, 
                           int surv = 1, int model = 1) 
  {
    return likelihood(vars[4],vars[0],vars[1],vars[2],vars[3],
                      est_norm,vflag,rescale,surv,model);
  }

  // ========================================================================

  void addRateShift(double t);

//  protected:
  // ========================================================================

  vector<double> times;     /* event times */
  vector<int> ttypes;       /* branching event type */
  int extant;               /* extant branches at t=0 */
  int maxExtant;            /* maximum number of extant branches */
  int nroot;                /* number of lineages at the root */
  string fn;                /* tree file name */
  int is_rev;               /* tree in forward or backwards format */
  // TreeNodeList nodes;    /* TreeNodes for extra information */
  int numBranches;
  vector< vector<int> > ids;
};

#endif
