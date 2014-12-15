#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
using namespace std;

#include "expoTree.h"
#include "Forest.h"
// #include "PSOLikelihood.h"
#include "TreeSwarm.h"
#include "pso/Swarm.h"
#include "pso/Point.h"

#ifdef INF_N
#define NUM_PARAMS 4
#else
#define NUM_PARAMS 5
#endif

extern "C" {
  void printHelp();
}

int parseArgs(int argc, char** argv);

// ===========================================================================

int Nmax = 1000;
double betaLo = 0.0;       /* infection rate */
double betaHi = 1.0;
double muLo = 0.0;         /* recovery rate */
double muHi = 0.1;
double psiLo = 0.4;        /* sampling rate */
double psiHi = 1.0;
double rhoLo = 0.0;        /* initial sampling probability */
double rhoHi = 1.0;
int root = 0;              /* add a root to the tree */
int SImodel = 1;           /* use non-saturating model */
int vflag = 0;             /* vebose flag */
int rescale = 1;           /* rescale probabilities in expomv */
int extant = 0;            /* number of extant species */
int maxExtant = 0;         /* max number of extant species */
int maxEvals = 1000;       /* max number of function evaluations */
int swarmSize = 20;        /* swarm size */
int optimAlg = 1;          /* lock variables */
double fixedRatio = -1.0;  /* fix mu-psi ratio */
long rngSeed = time(NULL); /* random number seed */

string fn = "";            /* input filename */
string out_fn = "";        /* output filename */
string hist_fn = "";       /* history filename */

// ===========================================================================

int main(int argc, char* argv[]) {
  if (parseArgs(argc,argv) == -1) return 0;
  if (optind == argc) {
    printHelp();
    return 0;
  }

  // =========================================================================

  ostream* oout = NULL;
  ostream* hist = NULL;

  if (out_fn != "") {
    if (out_fn == "-") oout = &cout;
    else oout = new ofstream(out_fn.c_str());
  }
  if (hist_fn != "") hist = new ofstream(hist_fn.c_str());

  // =========================================================================

  Forest trees(argc-optind,argv+optind);

  PSO::Parameters params(NUM_PARAMS);

  /*
#ifndef INF_N
  params.evalFunc = &likelihood;
#else
  params.evalFunc = &likelihood_inf;
#endif
  */

  params.type[0] = PSO::REAL;
  params.type[1] = PSO::REAL;
  params.type[2] = PSO::REAL;
  params.type[3] = PSO::REAL;

  params.lb[0] = betaLo; 
  params.lb[1] = muLo;
  params.lb[2] = psiLo;
  params.lb[3] = rhoLo;

  params.ub[0] = betaHi;
  params.ub[1] = muHi;
  params.ub[2] = psiHi;
  params.ub[3] = rhoHi;

  if ((optimAlg/8)  % 2 == 1) params.lockVar(0,betaLo);
  if ((optimAlg/4)  % 2 == 1) params.lockVar(1,muLo);
  if ((optimAlg/2)  % 2 == 1) params.lockVar(2,psiLo);
  if ((optimAlg)    % 2 == 1) params.lockVar(3,rhoLo);

  if (fixedRatio >= 0.0) params.lockVar(1,-fixedRatio);

#ifndef INF_N
  params.type[4] = PSO::INTEGER;
  params.lb[4] = (trees.maxExtant+1.);
  params.ub[4] = 1.*Nmax;
  if ((optimAlg/16) % 2 == 1) params.lockVar(4,Nmax);
#endif

  // custom evaluator creation:
  PSOParams pars;
  pars.SImodel = SImodel;
  pars.vflag = vflag;
  pars.rescale = rescale;
  pars.extant = extant;
  pars.maxExtant = maxExtant;
  pars.trees = &trees;
  pars.lockedRatio = fixedRatio;
  pars.tree = -1;
  params.evalParams = &pars;

  PSO::Point phi_p(NUM_PARAMS,1.5);
  PSO::Point phi_g(NUM_PARAMS,1.5);
  PSO::Point omega(NUM_PARAMS,0.7);
  PSO::TreeSwarm s(swarmSize,NUM_PARAMS,&params);

  s.seed_rng(rngSeed);
  s.setVars(phi_p,phi_g,omega);

  s.begin(argc,argv);

  if (s.is_master() && oout != NULL) {
    *oout << "# =====================\n"
          << "# Optimizing likelihood\n"
          << "# =====================\n";
#ifndef INF_N
    *oout << "# N        = [" << params.lb[4] << "," << params.ub[4] << "]\n";
#endif
    *oout << "# beta     = [" << betaLo << "," << betaHi << "]\n"
          << "# mu       = [" << muLo << "," << muHi << "]\n"
          << "# psi      = [" << psiLo << "," << psiHi << "]\n"
          << "# rho      = [" << rhoLo << "," << rhoHi << "]\n"
          << "# root     = " << root << "\n"
          << "# SImodel  = " << SImodel << "\n"
          << "# rescale  = " << rescale << "\n"
          << "# maxEvals = " << maxEvals << "\n"
          << "# lockPar  = " << optimAlg << "\n"
          << "# fixRatio = " << fixedRatio << "\n"
          << "# outfile  = " << out_fn << "\n" 
          << "# histfile = " << hist_fn << "\n"
          << "# times files:\n";
    for (size_t i(0); i < trees.size(); ++i) 
      *oout << "#   " << trees.at(i)->fn << "\n";
    *oout << "# \n"
          << "# extant species at time zero = " << extant << "\n"
          << "# max extant species = " << trees.maxExtant << "\n" << endl;
  }

  int cnt(0);
  while (1) {
    s.initialize(1);
    s.evaluate(vflag,NULL,hist);
#ifndef USE_MPI
    s.randomizeInf();
    if (s.best().fitness() == -INFINITY && cnt++ < 10) {
      cerr << "all inf! restarting..." << endl;
      for (int i(0); i < s.size(); ++i) cerr << s.at(i) << endl;
    } else {
      break;
    }
#else
    break;
#endif
  }

  // if (hist != NULL) s.display(hist);

#ifdef USE_MPI
  s.run_mpi(maxEvals,vflag,oout,hist);
#else
  s.run(maxEvals,0,vflag,oout,hist);
#endif

  if (s.is_master() && oout == NULL) cout << s.best() << endl;
  s.end();

  if (oout != &cout) delete oout;
  if (hist != NULL) delete hist;

  return 0;
}

// ===========================================================================

void printHelp() {
  printf("optim_pso\n");
  printf("  Find maximum likelihood value for a phylogenetic tree from an epidemic\n\n");
  printf("usage: optim_pso [-rvh] [-N <N>] [-b <minBeta>] [-B <maxBeta>] [-u <minMu>]\n");
  printf("                 [-s <minPsi>] [-S <maxPsi>] [-o <minRho>] [-O <maxRho>]\n");
  printf("                 [-l <SImodel>] [-t <fixPars>] timesFile1.txt timesFile2.txt ...\n\n");
  printf("  N : maximum population size\n");
  printf("  b : minimum infection rate 'beta'\n");
  printf("  B : maximum infection rate 'beta'\n");
  printf("  u : minimum recovery rate 'mu'\n");
  printf("  U : maximum recovery rate 'mu'\n");
  printf("  s : minimum sampling rate 'psi'\n");
  printf("  S : maximum sampling rate 'psi'\n");
  printf("  o : minimum initial sampling rate 'rho'\n");
  printf("  O : maximum initial sampling rate 'rho'\n");
  printf("  r : rescale staring vector after each iteration\n");
  printf("  l : model to use (0 = density independent; 1 = density dependent)\n");
  printf("  t : fixed some parameters\n");
  printf("  R : fix the ratio 'psi/(mu+psi)'\n");
  printf("  w : number of particles in the swarm\n");
  printf("  e : maximum number of evaluations\n");
  printf("  C : main output file (defaults to STDOUT)\n");
  printf("  H : history file\n");
  printf("  v : verbose (can be used multiple times)\n");
  printf("  h : print this help message\n");
}

// ===========================================================================

int parseArgs(int argc, char** argv) {
  int c;
  opterr = 0;
  while (1) {
    static struct option long_options[] = {
      {"verbose",     no_argument,       0, 'v'},
      {"maxPop",      required_argument, 0, 'N'},
      {"model",       required_argument, 0, 'l'},
      {"betaLo",      required_argument, 0, 'b'},
      {"betaHi",      required_argument, 0, 'B'},
      {"muLo",        required_argument, 0, 'u'},
      {"muHi",        required_argument, 0, 'U'},
      {"psiLo",       required_argument, 0, 's'},
      {"psiHi",       required_argument, 0, 'S'},
      {"rhoLo",       required_argument, 0, 'o'},
      {"rhoHi",       required_argument, 0, 'O'},
      {"maxEval",     required_argument, 0, 'e'},
      {"lockPar",     required_argument, 0, 't'},
      {"rescale",     required_argument, 0, 'r'},
      {"help",        no_argument,       0, 'h'},
      {"output",      required_argument, 0, 'C'},
      {"saveHist",    required_argument, 0, 'H'},
      {"swarmSize",   required_argument, 0, 'w'},
      {"lockRatio",   required_argument, 0, 'R'},
      {"rngSeed",     required_argument, 0, 'X'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc,argv,"N:l:vb:B:u:U:s:S:o:O:R:e:t:hC:H:w:X:r:",
        long_options,&option_index);
    if (c == -1) break;
    switch (c) {
      case 'N': Nmax = atoi(optarg); break;
      case 'l': SImodel = atoi(optarg); break;
      case 'v': ++vflag; break;
      case 'b': betaLo = atof(optarg); break;
      case 'B': betaHi = atof(optarg); break;
      case 'u': muLo = atof(optarg); break;
      case 'U': muHi = atof(optarg); break;
      case 's': psiLo = atof(optarg); break;
      case 'S': psiHi = atof(optarg); break;
      case 'o': rhoLo = atof(optarg); break;
      case 'O': rhoHi = atof(optarg); break;
      case 'R': fixedRatio = atof(optarg); break;
      case 'e': maxEvals = atoi(optarg); break;
      case 'r': rescale = atoi(optarg); break;
      case 't': optimAlg = atoi(optarg); break;
      case 'w': swarmSize = atoi(optarg); break;
      case 'C': out_fn = optarg; break;
      case 'H': hist_fn = optarg; break;
      case 'X': rngSeed = atoi(optarg); break;
      case 'h': printHelp(); return -1;
      default: break;
    }
  }

  return 0;
}

