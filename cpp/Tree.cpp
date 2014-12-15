#include "Tree.h"

double Tree::likelihood(vector<double>& K, vector<double>& beta, 
                        vector<double>& mu, vector<double>& psi, 
                        double rho, int est_norm, int vflag, 
                        int rescale, int surv, int model) 
{
  double fx(-INFINITY);

  switch (model) {
    case 0:
      break;
    case 2:
      fx = expoTreeSIR(K.size(),K.data(),beta.data(),mu.data(),psi.data(),
                       rho,times,ttypes,extant+nroot,est_norm,
                       vflag,rescale,nroot);
      break;
    case 1:
    default:
      fx = expoTreeEval(K,beta,mu,psi,rho,times,ttypes,extant,est_norm,vflag,rescale);
      break;
  }

  if (surv) {
    fx -= survival(K,beta,mu,psi,rho,est_norm,vflag,rescale,nroot,model);
  }

  return fx;
}

// ==========================================================================

double Tree::survival(vector<double>& K, vector<double>& beta, 
    vector<double>& mu, vector<double>& psi, double rho,
    int est_norm, int vflag, int rescale, int nroot, int model_type)
{
  double fx(-INFINITY);

  vector<int> survTypes(ttypes.size(),99);
  for (size_t i(0); i < ttypes.size(); ++i) {
    switch (ttypes[i]) {
      case 3:
        survTypes[i] = 3;
        break;
      default:
        break;
    }
  }

  switch (model_type) {
    case 0:
      break;
    case 2:
      fx = expoTreeSIRSurvival(K.size(),K.data(),beta.data(),mu.data(),
                      psi.data(),rho,times,ttypes,extant+nroot,est_norm,
                      vflag,rescale,nroot);
      break;

    case 1:
    default:
      fx = expoTreeSurvival(K,beta,mu,psi,rho,times,survTypes,
                            extant,est_norm,vflag,rescale);
  }
  return fx;
}

// ========================================================================

void Tree::reverse() {
  double mt = maxTime();
  for (size_t i(0); i < times.size(); ++i) {
    times[i] = mt-times[i];
  }
  std::reverse(times.begin(),times.end());
  std::reverse(ttypes.begin(),ttypes.end());
  std::reverse(ids.begin(),ids.end());
  is_rev = ! is_rev;
}

// ========================================================================

void Tree::addRateShift(double t) {
  // find position
  size_t pos(0);
  do {
    if (times[pos] > t) break;
    ++pos;
  } while (pos < times.size());
  vector<double>::iterator p1(times.begin());
  vector<int>::iterator p2(ttypes.begin());
  p1 += pos;
  p2 += pos;
  times.insert(p1,t);
  ttypes.insert(p2,20);
}

// ========================================================================

void Tree::readFromFile(string fn) {
  ifstream in(fn.c_str());
  if (! in.is_open()) {
    cerr << "Problem opening file '" << fn << "'. Aborting." << endl;
    return;
  }
  string input;
  double x1;
  int    x2;
  int    A;
  while (! getline(in,input).eof()) {
    if (input.length() == 0) continue;
    if (input[0] == '#') continue;
    istringstream istr(input);
    if (istr >> x1 >> x2) {
      times.push_back(x1);
      ttypes.push_back(x2);
      ids.push_back(vector<int>());
    } else {
      cerr << "Error reading times/type:" << endl << " > " << in << endl;
    }
    while (istr >> A) ids.back().push_back(A);
  }
  in.close();
  // make sure times are sorted
  vector<size_t> p(times.size());
  gsl_sort_index(p.data(),times.data(),1,times.size());
  // apply sort
  // gsl_permute(p.data(),times.data(),1,times.size());
  // gsl_permute_int(p.data(),ttypes.data(),1,ttypes.size());
  extant = 0;
  maxExtant = 0;
  numBranches = -1;
  int n(times.size());
  for (int i(n-1); i >= 0; --i) {
    switch (ttypes[i]) {
      case 1:
        ++extant;
        numBranches += 2;
        break;
      case 0:
      case 2:
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

