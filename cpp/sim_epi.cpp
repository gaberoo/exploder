#include <iostream>
#include <cmath>
using namespace std;

#include <mygsl/GSLRng.h>

int N;
double beta;
double mu;
double psi;
double rho;

int main(int argc, char* argv[]) {
  myGSL::Rng rng;

  N = atoi(argv[1]);
  beta = atof(argv[2]);
  mu = atof(argv[3]);
  psi = atof(argv[4]);
  rho = atof(argv[5]);

  int inf = rng.uniform_int(N-1)+1;
  int samp = 0;

  // Gillespie NRM
  double t = 0.0;
  double dt = 0.0;
  double total_rate = 0.0;
  double r;

  double rate[4] = { 0.0, 0.0, 0.0, 0.0 };

  while (inf > 0) {
    rate[0] = (inf-samp)*beta*(1-1.*inf/N);
    rate[1] = 2*samp*beta*(1-1.*inf/N);
    rate[2] = (inf-samp)*psi;
    rate[3] = (inf-samp)*mu;
    total_rate = rate[0] + rate[1] + rate[2] + rate[3];
    dt = -log(rng.uniform())/total_rate;
    t -= dt;
    r = rng.uniform();
    r -= rate[0]/total_rate;
    if (r < 0.0) {
      // infection event in one of the non-tree branches
      --inf;
    } else {
      r -= rate[1]/total_rate;
      if (r < 0.0) {
        // infection event in one of the tree branches
        --samp;
        --inf;
      } else {
        r -= rate[2]/total_rate;
        if (r < 0.0) {
          // sampling event
          ++samp;
        } else {
          // recovery in one of the non tree branches
          ++inf;
        }
      }
    }
    cout << t << " " << inf << " " << samp << endl;
  }

  return 0;
}
