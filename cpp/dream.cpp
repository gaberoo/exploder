// ===========================================================================
//
// DREAM sampler for estimating epidemic parameters
// ------------------------------------------------
//
// Gabriel E Leventhal
// Institute of Integrative Biology
// ETH Zurich
// Universitätstrasse 16
// 8092 Zürich
// Switzerland
//
// gabriel.leventhal@env.ethz.ch
// http://www.leventhal.ch
//
//
// DREAM algorithm:
//
// Vrugt, J. A., ter Braak, C. J. F., Diks, C. G. H., Robinson, B. A., Hyman, 
// J. M., Higdon, D., 2009. Accelerating Markov chain Monte Carlo simulation 
// by differential evolution with self-adaptive randomized subspace sampling. 
// International Journal of Nonlinear Sciences and Numerical Simulation 
// 10 (3), 273-290. DOI: 10.1515/IJNSNS.2009.10.3.273
//
// ===========================================================================

#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <streambuf>
#include <vector>
#include <algorithm>
#include <map>
#include <getopt.h>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics.h>

#include <rapidjson/document.h>

#include "array.h"
#include "expoTree.h"
#include "Forest.h"
#include "gelman_rubin.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

#ifdef PF
#include <omp.h>
#include "mctraj/pfLik.h"
#endif

// ===========================================================================

void printHelp();
int parseArgs(int argc, char** argv);

#ifdef USE_MPI
double lik_master(Forest& trees, const vector<double>& vars,
                  int SImodel, int vflag, int rescale);
void lik_slave(Forest& trees, vector<double>& vars, 
               int SImodel, int vflag, int rescale, int survival);
#endif

// ===========================================================================

void check_outliers(int t, Array2D<double>& lik, vector<double>& meanlik,
                    vector<bool>& outliers) 
{
  int t0 = t/2;
  int numChains = lik.n_y();
  double Q1;
  double Q3;
  double IQR;
  double UR;
  vector<double> liksrt(numChains);
  if ((int) meanlik.size() < numChains) meanlik.resize(numChains,-INFINITY);
  for (int i(0); i < numChains; ++i) {
    meanlik[i] = gsl_stats_mean(lik.pt(t0,i),numChains,t-t0);
    liksrt[i] = meanlik[i];
  }
  sort(liksrt.begin(),liksrt.end());
  Q1 = gsl_stats_quantile_from_sorted_data(liksrt.data(),1,liksrt.size(),0.25);
  Q3 = gsl_stats_quantile_from_sorted_data(liksrt.data(),1,liksrt.size(),0.75);
  IQR = Q3 - Q1;
  UR = Q1 - 2*IQR;
  for (int i(0); i < numChains; ++i) {
    outliers[i] = false;
    if (meanlik[i] < UR) outliers[i] = true;
  }
}

// ===========================================================================

void gen_CR(const gsl_rng* rng, const vector<double>& pCR, 
    Array2D<int>& CRm, vector<unsigned>& L) 
{
  size_t numChains(CRm.n_x());
  size_t loopSteps(CRm.n_y());
  size_t size(numChains*loopSteps);
  if (L.size() < pCR.size()) {
    cerr << "Bad array size." << endl;
    return;
  }
  // pick candidates for each crossover value
  gsl_ran_multinomial(rng,pCR.size(),size,pCR.data(),L.data());
  size_t j(0);
  for (size_t m(0); m < L.size(); ++m) {
    for (size_t i(0); i < L[m]; ++i) {
      CRm[j] = m+1;
      ++j;
    }
  }
  gsl_ran_shuffle(rng,CRm.pt(0,0),size,sizeof(int));
}

// ===========================================================================

// ===========================================================================

double Nmin = 0;           /* minimum total population size */
double Nmax = 1000;        /* maximum total population size */
int root = 0;              /* add a root to the tree */
int SImodel = 1;           /* use non-saturating model */
int n = 0;                 /* number of times */
int vflag = 0;             /* vebose flag */
int rescale = 1;           /* rescale probabilities in expomv */
int extant = 0;            /* number of extant species */
int maxExtant = 0;         /* max number of extant species */
int maxEvals = 100000;     /* max number of function evaluations */
int optimAlg = 1;          /* lock variables */
double fixedRatio = -1.0;  /* fix mu-phi ratio */
char jointLikelihood = 'j';   /* forest likelihood type */
int numChains = 5;
int survival = 1;          /* correct for survival probability of the tree */
int modelType = 1;
string json_input = "[]";

string fn = "";            /* input filename */
string out_fn = "";        /* output filename */
int appendFile = 0;        /* continue from previous state */
int report_interval = 1;   /* report interval for state */
int diagnostics = 0;       /* report diagnostics at the end of the run */
int reverseTrees = 1;
int saveTraj = 0;

int burnIn = 0;            /* number of steps for which to run an adaptive proposal size */

// DREAM variables
double noise = 0.05;
double bstar_zero = 1e-3;
int collapseOutliers = 1;
int gelmanEvals = 100;
int loopSteps = 10;
double scaleReductionCrit = 1.01;
int numPars = 5;
int deltaMax = 2;
int pCR_update = 1;
int nCR = 3;
double reenterBurnin = 0.2;

// MPI variables
int mpi_rank = 0;
int mpi_ntasks = 0;
clock_t calc_time = 0;
clock_t run_time = 0;
clock_t last_time = 0;

#ifdef PF
// Particle Filtering variables
int num_particles = 1000;
vector<ostream*> traj_out;
vector<ostream*> branch_out;
double pfLikelihood(Tree& tree, const vector<double>& vars,
                    int num_particles, gsl_rng** rng, Trajectory* ttraj = NULL) 
{
  MCTraj::SISModel::EpiPars pars = { vars[4], vars[0], vars[1], vars[2], vars[3] };
  SEISModel::EpiPars seis_pars;

  EpiState* es = NULL;
  Model* mpt = NULL;

  switch (modelType) {
    case 1:
    default:
      mpt = new SIS(&pars);
      es = new EpiState(SISModel::nstates);
      (*es)[0] = (int) pars.N;
      break;
    case 2:
      mpt = new SIR(&pars);
      es = new EpiState(SISModel::nstates);
      (*es)[0] = (int) pars.N;
      break;
    case 3:
      seis_pars.N = vars[4];
      seis_pars.beta = vars[0];
      seis_pars.mu = vars[1];
      seis_pars.psi = vars[2];
      seis_pars.rho = vars[3];
      seis_pars.gamma = vars[5];
      mpt = new SEIS(&seis_pars);
      es = new EpiState(SEISModel::nstates);
      (*es)[0] = (int) seis_pars.N;
      break;
  }

  (*es)[0] = (int) pars.N;

  double lik = -INFINITY;
  if (! tree.is_rev) tree.reverse();
  es->init_branches(tree.max_id()+1);

  lik = MCTraj::pfLik(mpt,*es,tree,num_particles,rng,vflag,ttraj);

  delete mpt;

  return lik;
}
#endif

// ===========================================================================

int main(int argc, char* argv[]) {
#ifdef PF
  // omp_set_num_threads(1);
  // collapseOutliers = 0;
#endif

#ifdef USE_MPI
  MPI::Init(argc,argv);
  mpi_rank = MPI::COMM_WORLD.Get_rank();
  mpi_ntasks = MPI::COMM_WORLD.Get_size();
#endif

  // =========================================================================

  if (parseArgs(argc,argv) == -1) return 0;
  if (optind == argc) {
    cerr << "Please provide at least one tree file.\n" << endl;
    printHelp();
    return 0;
  }

  // =========================================================================
  // lock certain parameters

  // =========================================================================

  /* Reading priors from JSON */
  rapidjson::Document jpars;
  rapidjson::Value::Member* m1;
  rapidjson::SizeType _rj; 

  jpars.Parse<0>(json_input.c_str());
  assert(jpars.IsObject());

  rapidjson::Value::Member* _d = jpars.FindMember("vars");
  assert(_d != 0);

  rapidjson::Value& d = _d->value;

  int numPars = d.Size();
  vector<double> varInit(numPars,-1.0);
  vector<double> varLo(numPars,0.0);
  vector<double> varHi(numPars,INFINITY);
  vector<bool> lockVar(numPars,false);
  vector<string> varNames(numPars,"");
  map<string,int> nameMap;
  int realPars = numPars;

  // get names
  for (rapidjson::SizeType i = 0; i < d.Size(); ++i) {
    rapidjson::Value& a = d[i];
    m1 = a.FindMember("name"); 
    if (m1 != 0) {
      assert(m1->value.IsString());
      varNames[i] = m1->value.GetString();
      nameMap.insert(make_pair(varNames[i],i));
    }
  }

  for (rapidjson::SizeType i = 0; i < d.Size(); ++i) {
    rapidjson::Value& a = d[i];

    m1 = a.FindMember("limits"); 
    assert(m1 != 0); 
    assert(m1->value.IsArray());
    assert(m1->value.Size() == 2);
    _rj = 0; varLo[i] = m1->value[_rj].GetDouble();
    _rj = 1; varHi[i] = m1->value[_rj].GetDouble();

    m1 = a.FindMember("lock");
    if (m1 != 0) {
      assert(m1->value.IsBool());
      lockVar[i] = m1->value.GetBool();
      if (lockVar[i]) --realPars;
    }

    m1 = a.FindMember("fixed_to");
    if (m1 != 0) {

    }
  }

  m1 = jpars.FindMember("save_trajectory");
  if (m1 != 0) {
    assert(m1->value.IsBool());
    saveTraj = m1->value.GetBool();
  }

  m1 = jpars.FindMember("prefix");
  if (m1 != 0) {
    assert(m1->value.IsString());
    out_fn = m1->value.GetString();
  }

  m1 = jpars.FindMember("num_chains");
  if (m1 != 0) {
    assert(m1->value.IsInt());
    numChains = m1->value.GetInt();
  }

  m1 = jpars.FindMember("max_evals");
  if (m1 != 0) {
    assert(m1->value.IsInt());
    maxEvals = m1->value.GetInt();
  }

  m1 = jpars.FindMember("model_type");
  if (m1 != 0) {
    assert(m1->value.IsInt());
    modelType = m1->value.GetInt();
  }

  m1 = jpars.FindMember("burn_in");
  if (m1 != 0) {
    assert(m1->value.IsInt());
    burnIn = m1->value.GetInt();
  }

  m1 = jpars.FindMember("gelman_evals");
  if (m1 != 0) {
    assert(m1->value.IsInt());
    gelmanEvals = m1->value.GetInt();
  }

  m1 = jpars.FindMember("survival");
  if (m1 != 0) {
    assert(m1->value.IsInt());
    survival = m1->value.GetInt();
  }

  if (fixedRatio >= 0.0) lockVar[1] = true;
  deltaMax = (realPars-1)/2;

  // =========================================================================
  // load trees

  Forest trees(argc-optind,argv+optind);
#ifdef PF
  if (reverseTrees) {
    for (size_t ii(0); ii < trees.size(); ++ii) trees[ii]->reverse();
  }
#endif

  // =========================================================================
  // START OF ALGORITHM
  // =========================================================================
  
#ifdef PF
  vector<MCTraj::Trajectory> traj(numChains);
  vector<MCTraj::Trajectory> traj_proposal(numChains);
#endif
  
  int inBurnIn = (burnIn > 0);
  int burnInStart = 0;

  // set initial values

  for (int i(0); i < numPars; ++i) {
    if (varInit[i] < 0.0) varInit[i] = varLo[i];
  }
  if (fixedRatio >= 0.0) varInit[1] = varInit[2]*(1./fixedRatio-1.);

  int _nvar = nameMap["N"];
  if (! lockVar[_nvar]) {
    if (varLo[4] < trees.maxExtant) varLo[4] = trees.maxExtant;
  }

  cerr << "Beginning DREAM:" << endl;
  for (int i = 0; i < numPars; ++i) {
    cerr << setw(16) << varNames[i] << " = [" 
         << setw(6) << varLo[i] << "," 
         << setw(6) << varHi[i] << "]";
    if (lockVar[i]) cerr << " *";
    cerr << endl;
  }

  // int survival = (lockVar[2] && varLo[2] == 0.0 && lockVar[3] && varLo[3] == 0.0) ? 0 : 1;

  int genNumber(0);
  if (mpi_rank == 0) {
    // initialize random number generator
    gsl_rng* rng = NULL;
#ifdef PF
    gsl_rng** _rng = NULL;
    size_t r = 0;
    size_t nthreads = omp_get_max_threads();
    long int seed = time(NULL);
    _rng = (gsl_rng**) malloc(nthreads*sizeof(gsl_rng*));
    for (r = 0; r < nthreads; ++r) {
      _rng[r] = gsl_rng_alloc(gsl_rng_taus2);
      gsl_rng_set(_rng[r],seed+r);
    }
    rng = _rng[0];
#else
    rng = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng,time(NULL));
#endif

    // MCMC chains
    Array3D<double> state(maxEvals,numChains,numPars);
    Array2D<double> lik(maxEvals+1,numChains);

    vector< vector<double> > proposal(numChains,vector<double>(numPars,0.0));
    vector<double> scaleReduction(numPars,0.0);
    vector<double> pCR(nCR,1./nCR);

    // =========================================================================
    // read previous state

    int prevLines(0);
    if (appendFile && (out_fn != "" || out_fn != "-")) {
      cerr << "Restoring previous state... ";
      prevLines = maxEvals;
      for (int i(0); i < numChains; ++i) {
        cerr << i << " ";
        int line(0);
        ostringstream chain_fn("");
        chain_fn << out_fn << "." << i << ".txt";
        ifstream ifile(chain_fn.str().c_str());
        if (! ifile) {
          appendFile = 0;
          prevLines = 0;
          cerr << "files don't exists: " << chain_fn.str() << endl;
          break;
        } else {
          string input;
          while (! getline(ifile,input).eof()) {
            if (line >= maxEvals) break;
            if (input.length() == 0) continue;
            if (input[0] == '#') continue;
            istringstream istr(input);
            for (int j = 0; j < numPars; ++j) istr >> state(line,i,j); 
            istr >> lik(line,i) >> inBurnIn >> genNumber;
            for (int j(0); j < nCR; ++j) istr >> pCR[j];
            ++line;
          }
          ifile.close();
          if (prevLines > line) prevLines = line-1;
        }
        for (int j(0); j < nCR; ++j) cerr << pCR[j] << " ";
        cerr << endl;
      }
    }

    // =========================================================================
    // open output file

    ostringstream chain_fn;
    vector<ostream*> oout(numChains,&cout);
    ios_base::openmode fmode = (appendFile) ? (ios_base::out | ios_base::app) : ios_base::out;

#ifdef PF
    if (saveTraj) {
      traj_out.resize(numChains);
      branch_out.resize(numChains);
    }
#endif

    if (out_fn != "" && out_fn != "-") {
      for (int i(0); i < numChains; ++i) {
        cerr << "opening output file " << i << "...";

        chain_fn.str("");
        chain_fn << out_fn << "." << i << ".txt";
        oout[i] = new ofstream(chain_fn.str().c_str(),fmode);
        oout[i]->setf(ios::scientific,ios::floatfield);
        oout[i]->precision(12);
        cerr << "done." << endl;
        if (appendFile) *oout[i] << "# --- Resuming DREAM ---" << endl;

#ifdef PF
        // =====================================================================
        // open PF files
        if (saveTraj) {
          chain_fn.str("");
          chain_fn << out_fn << "." << i << ".traj";
          traj_out[i] = new ofstream(chain_fn.str().c_str(),fmode);
          if (appendFile) *traj_out[i] << "# --- Resuming DREAM ---" << endl;

          chain_fn.str("");
          chain_fn << out_fn << "." << i << ".branch";
          branch_out[i] = new ofstream(chain_fn.str().c_str(),fmode);
          if (appendFile) *branch_out[i] << "# --- Resuming DREAM ---" << endl;
        }
#endif
      }
    }
 
    // =========================================================================
    // Initialize with latin hypbercube sampling if not a resumed run

    int do_calc(1);

    if (! appendFile) {
      cerr << "Latin hypercube sampling..." << endl;

      vector< vector<int> > samples(numPars,vector<int>(numChains,0));
      // Array2D<int> samples(numPars,numChains);
      for (int j(0); j < numPars; ++j) {
        for (int i(0); i < numChains; ++i) samples[j][i] = i;
        gsl_ran_shuffle(rng,samples[j].data(),numChains,sizeof(int));
      }

      double randPos;
      double intervalSize;
      for (int i(0); i < numChains; ++i) {
        for (int j(0); j < numPars; ++j) {
          if (lockVar[j]) randPos = varLo[j];
          else {
            intervalSize = (varHi[j]-varLo[j])/numChains;
            randPos = varLo[j] + samples[j][i]*(varHi[j]-varLo[j])/numChains
                        + intervalSize*gsl_rng_uniform(rng);
          }
          if (randPos < varLo[j]) randPos = varLo[j];
          else if (randPos > varHi[j]) randPos = varHi[j];
          state(0,i,j) = randPos;
        }
        if (fixedRatio >= 0.0) state(0,i,1) = state(0,i,2)*(1./fixedRatio-1.);
      }

      for (int i(0); i < numChains; ++i) {
        do_calc = 1;
        for (int j(0); j < numPars; ++j) {
          if (state(0,i,j) < varLo[j] || state(0,i,j) > varHi[j]) {
            do_calc = 0;
            break;
          }
        }
        if (do_calc) {
          vector<double> vars(numPars);
          for (int j(0); j < numPars; ++j) vars[j] = state(0,i,j);
#ifdef USE_MPI
          lik(0,i) = lik_master(trees,vars,SImodel,vflag,rescale);
#else
#ifdef PF
          lik(0,i) = 0.0;
          for (size_t ii(0); ii < trees.size(); ++ii) {
            double x = pfLikelihood(*trees[ii],vars,num_particles,_rng,&traj[i]);
            lik(0,i) += x;
          }
#else
          lik(0,i) = trees.jointLikelihood(vars,SImodel,vflag,rescale,survival);
#endif
#endif
          if (vflag) cout << "Chain " << i << " likelihood = " << lik(0,i) << endl;
        } else lik(0,i) = -INFINITY;
      }

      // replace chains with infinite likelihood with random samples
      for (int i(0); i < numChains; ++i) {
        while (lik(0,i) == -INFINITY || lik(0,i) != lik(0,i)) {
          cerr << "Likelihood of chain = INF. Resampling intial parameters..." << endl;
          // choose random parameters
          for (int j(0); j < numPars; ++j) {
            if (lockVar[j]) randPos = varLo[j];
            else {
              randPos = varLo[j] + (varHi[j]-varLo[j])*gsl_rng_uniform(rng);
            }
            if (randPos < varLo[j]) randPos = varLo[j];
            else if (randPos > varHi[j]) randPos = varHi[j];
            state(0,i,j) = randPos;
          }
          if (fixedRatio >= 0.0) state(0,i,1) = state(0,i,2)*(1./fixedRatio-1.);
          vector<double> vars(numPars);
          for (int j(0); j < numPars; ++j) vars[j] = state(0,i,j);
#ifdef USE_MPI
          lik(0,i) = lik_master(trees,vars,SImodel,vflag,rescale);
#else
#ifdef PF
          lik(0,i) = 0.0;
          for (size_t ii(0); ii < trees.size(); ++ii) {
            lik(0,i) += pfLikelihood(*trees[ii],vars,num_particles,_rng,&traj[i]);
          }
#else
          lik(0,i) = trees.jointLikelihood(vars,SImodel,vflag,rescale,survival);
#endif
#endif
          cerr << "New likelihood = " << lik(0,i) << endl;
        }
      }

      // save initial state of each chain
      for (int i(0); i < numChains; ++i) {
        for (int j = 0; j < numPars; ++j) *oout[i] << state(0,i,j) << " ";
        *oout[i]  << lik(0,i) << " " << inBurnIn << " " << 0 << " ";
        for (int j(0); j < nCR; ++j) *oout[i] << pCR[j] << " ";
        *oout[i] << endl;
      }
    }

    // =========================================================================
 
    double gamma;
    int r1, r2;
    double crossRate(1.0);
    int curRun(0);
    int gammaGeneration(0);
    int delta;
    int numAccepted(0);
    int ireport(0);
    vector<int> acceptStep(numChains,0);

    vector<double> pairDiff(numPars);
    vector<double> e(numPars);
    vector<double> epsilon(numPars);
    vector<bool>   updatePar(numPars,false);
    vector<int>    updateDim(numChains,realPars);
    vector<double> step(numPars);

    vector<double> bstar(numPars,bstar_zero);
    bstar[4] = 1;

    // ======================================================================
    // RUN MCMC
    // ======================================================================
    
    vector<unsigned> L(nCR,0);      // candidates for crossover
    vector<unsigned> totalSteps(nCR,0);
    Array2D<int> CRm(numChains,loopSteps);
    vector<double> delta_tot(nCR,1.0);
    vector<double> sd(numPars,0.0);
    vector<double> delta_normX(numChains,0.0);
    double delta_sum(0.0);
    double pCR_sum(0.0);

    for (int t(prevLines+1); t < maxEvals; ++t) {
      // beginning of loop, generate crossover probabilities
      if (genNumber == 0) { 
        gen_CR(rng,pCR,CRm,L); 
        numAccepted = 0;
      }

      for (int i(0); i < numChains; ++i) {
        if (deltaMax > 1) delta = gsl_rng_uniform_int(rng,deltaMax)+1;
        else delta = 1;

        // generate proposal
        for (int j(0); j < numPars; ++j) proposal[i][j] = state(t-1,i,j);
        if (gammaGeneration++ == 5) {
          for (int j(0); j < numPars; ++j) {
            if (lockVar[j]) continue;
            gammaGeneration = 0;
            do {
              r1 = gsl_rng_uniform_int(rng,numChains-1);
              r2 = gsl_rng_uniform_int(rng,numChains-1);
              if (r1 >= i) ++r1;
              if (r2 >= i) ++r2;
            } while (r1 == r2 && numChains > 2);
            step[j] = state(t-1,r1,j) - state(t-1,r2,j);
            proposal[i][j] = state(t-1,i,j) + step[j];
          }
          // (TODO) if step == 0 ==> use Cholskey decomposition!!
        } else {
          // pick pairs
          vector<int> r1(delta,0);
          vector<int> r2(delta,0);
          for (int a(0); a < delta; ++a) {
            r1[a] = gsl_rng_uniform_int(rng,numChains-1);
            r2[a] = gsl_rng_uniform_int(rng,numChains-1);
            if (r1[a] >= i) ++r1[a];
            if (numChains > 2) {
              while (r2[a] == r1[a] || r2[a] == i) ++r2[a] %= numChains;
            }
          }
          updateDim[i] = realPars;
          for (int j(0); j < numPars; ++j) {
            if (! lockVar[j]) {
              updatePar[j] = true;
              // calculate random values
              e[j] = noise*(2.*gsl_rng_uniform(rng)-1.);
              epsilon[j] = gsl_ran_gaussian(rng,bstar[j]);
              // compute pairwise comparisons
              pairDiff[j] = 0.0;
              for (int a(0); a < delta; ++a) {
                if (r1[a] != r2[a]) pairDiff[j] += state(t-1,r1[a],j) - state(t-1,r2[a],j);
              }
              // check for crossover events
              crossRate = 1.*CRm(i,genNumber)/nCR;
              if (crossRate < 1.0) {
                if (gsl_rng_uniform(rng) < 1.0-crossRate) {
                  step[j] = 0.0;
                  updatePar[j] = false;
                  --updateDim[i];
                }
              }
            }
          }
          if (updateDim[i] > 0) {
            gamma = 2.38/sqrt(2.0*updateDim[i]*delta)/sqrt(trees.size());
            for (int j(0); j < numPars; ++j) {
              if (updatePar[j]) {
                // calculate step for this dimension
                step[j] = (1+e[j])*gamma*pairDiff[j] + epsilon[j];
                // update proposal
                proposal[i][j] = state(t-1,i,j) + step[j];
              } else {
                proposal[i][j] = state(t-1,i,j);
              }
            }
          } else {
            for (int j(0); j < numPars; ++j) proposal[i][j] = state(t-1,i,j);
          }
        }
      }

      // calculate likelihood and acceptance probability
      for (int i(0); i < numChains; ++i) {
        if (fixedRatio >= 0.0) proposal[i][1] = proposal[i][2]*(1./fixedRatio-1.);
      }

      for (int i(0); i < numChains; ++i) {
        if (updateDim[i] > 0) {
          do_calc = 1;
          for (int j(0); j < numPars; ++j) {
            if (!lockVar[j]) {
              if (proposal[i][j] < varLo[j] || proposal[i][j] > varHi[j]) {
                do_calc = 0;
                break;
              }
            }
          }
          if (do_calc) {
            if (vflag) {
              cout << "Chain " << i << ": ";
              cout << "N=" << proposal[i][4] << ", beta=" << proposal[i][0] << ", ";
              cout << "mu=" << proposal[i][1] << ", psi=" << proposal[i][2] << ", ";
              cout << "psi=" << proposal[i][3] << flush;
            }
#ifdef USE_MPI
            lik(t,i) = lik_master(trees,proposal[i],SImodel,vflag,rescale);
#else
#ifdef PF
            lik(t,i) = 0.0;
            for (size_t ii(0); ii < trees.size(); ++ii) {
              lik(t,i) += pfLikelihood(*trees[ii],proposal[i],num_particles,_rng,&traj_proposal[i]);
            }
            while (lik(t,i) == 0.0) {
              cerr << "Likelihood is zero!!" << endl;
              lik(t,i) = 0.0;
              for (size_t ii(0); ii < trees.size(); ++ii) {
                lik(t,i) += pfLikelihood(*trees[ii],proposal[i],num_particles,_rng,&traj_proposal[i]);
              }
            }
#else
            lik(t,i) = trees.jointLikelihood(proposal[i],SImodel,vflag,rescale,survival);
#endif
#endif
            if (vflag) cout << ". Likelihood = " << lik(t,i) << endl;
          } else lik(t,i) = -INFINITY;
        } else {
          for (int j(0); j < numPars; ++j) proposal[i][j] = state(t-1,i,j);
          lik(t,i) = lik(t-1,i);
        }
      }

      for (int i(0); i < numChains; ++i) {
        if (lik(t,i) >= lik(t-1,i)) acceptStep[i] = 1;
        else if (lik(t,i) == -INFINITY) acceptStep[i] = 0;
        else {
          if (log(gsl_rng_uniform(rng)) < lik(t,i)-lik(t-1,i)) acceptStep[i] = 1;
          else acceptStep[i] = 0;
        }
        if (acceptStep[i]) {
          ++numAccepted;
          for (int j(0); j < numPars; ++j) state(t,i,j) = proposal[i][j];
#ifdef PF
          traj[i] = traj_proposal[i];
#endif
        } else {
          for (int j(0); j < numPars; ++j) state(t,i,j) = state(t-1,i,j);
          lik(t,i) = lik(t-1,i);
        }
      }

      // update pCR if in burn-in phase
      if (inBurnIn && pCR_update) {
        // get standard deviations between the chains
        for (int j(0); j < numPars; ++j) {
          sd[j] = gsl_stats_sd(state.pt(t,0,j),numPars,numChains);
//          if (! lockVar[j] && sd[j] == 0.0) {
//            bstar[j] *= 2;
//            cerr << "Variable " << j << " has collapsed. Increasing stochasticity to " << bstar[j] << "." << endl;
//          }
        }
        for (int i(0); i < numChains; ++i) {
          if (acceptStep[i]) {
            // get Euclidian distance
            delta_normX[i] = 0.0;
            for (int j(0); j < numPars; ++j) {
              if (! lockVar[j] && sd[j] > 0.0) {
                delta_normX[i] += gsl_pow_2((state(t,i,j)-state(t-1,i,j))/sd[j]);
              }
            }
            delta_tot[CRm(i,genNumber)-1] += delta_normX[i];
          }
        }
      }

      if (++genNumber >= loopSteps) {
        genNumber = 0;

        if (inBurnIn && pCR_update) {
          // get total delta
          delta_sum = 0.0;
          for (int m(0); m < nCR; ++m) {
            delta_sum += delta_tot[m];
            totalSteps[m] += L[m];
          }
          if (delta_sum > 0.0) {
            pCR_sum = 0.0;
            for (int m(0); m < nCR; ++m) {
              pCR[m] = (delta_tot[m]/delta_sum) / totalSteps[m];  // relative jump size per step
              pCR_sum += pCR[m];
            }
            for (int m(0); m < nCR; ++m) {
              pCR[m] /= pCR_sum;
            }
          }
        }

        if (collapseOutliers && t < reenterBurnin*maxEvals) {
          // remove outlier chains
          vector<double> meanlik(numChains,-INFINITY);
          vector<bool> outliers(numChains,false);
          check_outliers(t,lik,meanlik,outliers);
          int best_chain = gsl_stats_max_index(meanlik.data(),1,numChains);
          for (int i(0); i < numChains; ++i) {
            if (outliers[i] && i != best_chain) {
              // chain is an outlier
              lik(t,i) = lik(t,best_chain);
              for (int j(0); j < numPars; ++j) state(t,i,j) = state(t,best_chain,j);
              if (! inBurnIn && burnIn > 0) {
                cerr << "[" << t << "] Outlier chain detected [" << i << "] outside of burn in."
                     << " Moving to chain " << best_chain << " and re-entering burn in." << endl;
                burnInStart = t;
                curRun = 0;
                inBurnIn = 1;
              }
            }
          }
        }

        // check if in burn in period
        if (! inBurnIn) {
          // calculate Gelman-Rubin convergence
          if (curRun++ >= gelmanEvals) {
            cout << "[" << t << "] performing convergence diagnostics:";
            gelman_rubin(state,scaleReduction,lockVar,burnInStart+burnIn,t);
            // estimate variance
            int exitLoop(numPars);
            cout << " GR (it " << t << "): ";
            for (int j(0); j < numPars; ++j) {
              if (! lockVar[j]) {
                // scale reduction factor
                if (scaleReduction[j] < scaleReductionCrit) --exitLoop;
                cout << "[" << j << "]" << scaleReduction[j] << " ";
              } else --exitLoop;
            }
            cout << endl;
            // check for convergence
            if (exitLoop <= 0) break;
            // reset counter
            curRun = 0;
          }
        } else {
          // check if burn-in is finished
          if (t >= burnInStart + burnIn) {
            inBurnIn = 0;
            cerr << "[" << t << "] exiting burn in." << endl;
          }
        }
      }

      ++ireport;
      if (ireport >= report_interval) {
        ireport = 0;
        for (int i(0); i < numChains; ++i) {
          for (int j = 0; j < numPars; ++j) *oout[i] << state(t,i,j) << " ";
          *oout[i] << lik(t,i) << " " << (t < burnInStart+burnIn) << " " << genNumber << " ";
          for (int j(0); j < nCR; ++j) *oout[i] << pCR[j] << " ";
          *oout[i] << endl;

#ifdef PF
          if (saveTraj) {
            traj[i].printBranches(*branch_out[i]) << endl;
            traj[i].printFromFirst(*traj_out[i]) << endl;
          }
#endif
        }
      }
    }

#ifdef USE_MPI
    for (int k = 1; k < mpi_ntasks; ++k) {
      if (vflag) fprintf(stderr,"Sending done signal to process %d.\n",k);
      MPI::COMM_WORLD.Send(0,0,MPI::INT,k,0);
    }
#endif

#ifdef PF
    for (r = 0; r < nthreads; ++r) {
      gsl_rng_free(_rng[r]);
      _rng[r] = NULL;
    }
    free(_rng);
    rng = NULL;
#else
    gsl_rng_free(rng);
    rng = NULL;
#endif

    for (int i(0); i < numChains; ++i) delete oout[i];
#ifdef PF
    if (saveTraj) {
      for (int i(0); i < numChains; ++i) {
        delete traj_out[i];
        delete branch_out[i];
      }
    }
#endif
  } else {
#ifdef USE_MPI
    vector<double> newVar(numPars,0.0);
    lik_slave(trees,newVar,SImodel,0,rescale,survival);
#endif
  }

#ifdef USE_MPI
  last_time = clock();
  fprintf(stdout,"Calcluation efficiency (Slave %d): %f\n",
          mpi_rank,1.*calc_time/last_time);
  MPI::Finalize();
#endif

  return EXIT_SUCCESS;
}

// ===========================================================================

void printHelp() {
  printf("mcmc: Markov chain Monte Carlo for density dependent trees\n\n");
  printf("usage: mcmc -N <population_size> -b <beta> -u <mu> -s <psi>\n");
  printf("               -l <SImodel> -f <times_file>\n");
  printf("               [-rvh]\n\n");
  printf("  N : starting population size\n");
  printf("  b : starting infection rate 'beta'\n");
  printf("  u : starting recovery rate 'mu'\n");
  printf("  r : rescale staring vector after each iteration\n");
  printf("  s : starting sampling rate 'psi'\n");
  printf("  l : model to use (0 = density independent; 1 = density dependent)\n");
  printf("  f : file with time samples\n");
  printf("  t : fixed parameters\n");
  printf("  v : verbose (can be used multiple times)\n");
  printf("  h : print this help message\n");
}

// ===========================================================================

int parseArgs(int argc, char** argv) {
  int c;
  opterr = 0;
  ifstream in;
  while (1) {
    static struct option long_options[] = {
      {"append",        no_argument,       0, 'A'},
      {"collapseOut",   required_argument, 0, 'c'},
      {"output",        required_argument, 0, 'C'},
      {"diagnostics",   no_argument,       0, 'D'},
      {"maxEval",       required_argument, 0, 'e'},
      {"reenterBurnin", required_argument, 0, 'E'},
      {"gelmanEvals",   required_argument, 0, 'G'},
      {"gelmanCrit",    required_argument, 0, 'g'},
      {"help",          no_argument,       0, 'h'},
      {"report",        required_argument, 0, 'I'},
      {"jointType",     required_argument, 0, 'J'},
      {"model",         required_argument, 0, 'l'},
      {"nparticles",    required_argument, 0, 'm'},
      {"maxPop",        required_argument, 0, 'M'},
      {"fixedRatio",    required_argument, 0, 'R'},
      {"rescale",       required_argument, 0, 'r'},
      {"psiLo",         required_argument, 0, 's'},
      {"psiHi",         required_argument, 0, 'S'},
      {"lockPar",       required_argument, 0, 't'},
      {"verbose",       no_argument,       0, 'v'},
      {"burnIn",        required_argument, 0, 'X'},
      {"json",          required_argument, 0, 'z'},
      {"survival",      required_argument, 0, 'x'},
      {0, 0, 0, 0}
    };
    int option_index = 0;
    c = getopt_long(argc,argv,"a:Ac:C:De:E:G:g:hI:J:l:m:M:R:r:t:T:vWw:x:X:z:",
        long_options,&option_index);
    if (c == -1) break;
    switch (c) {
      case 'l': SImodel = atoi(optarg); break;
      case 'v': ++vflag; break;
      case 'R': fixedRatio = atof(optarg); break;
      case 'e': maxEvals = atoi(optarg); break;
      case 'E': reenterBurnin = atof(optarg); break;
      case 'r': rescale = atoi(optarg); break;
      case 'T': modelType = atoi(optarg); break;
      case 't': optimAlg = atoi(optarg); break;
      case 'C': out_fn = optarg; break;
      case 'h': printHelp(); return -1;
      case 'I': report_interval = atoi(optarg); break;
      case 'D': ++diagnostics; break;
      case 'X': burnIn = atoi(optarg); break;
      case 'J': jointLikelihood = optarg[0]; break;
      case 'A': appendFile = 1; break;
      case 'G': gelmanEvals = atoi(optarg); break;
      case 'g': scaleReductionCrit = GSL_MAX(atof(optarg),1.0); break;
      case 'c': collapseOutliers = atoi(optarg); break;
      case 'x': survival = atoi(optarg); break;
      case 'z':
        in.open(optarg);
        in.seekg(0,ios::end);
        json_input.reserve(in.tellg());
        in.seekg(0,ios::beg);
        json_input.assign(istreambuf_iterator<char>(in), 
                          istreambuf_iterator<char>());
        break;
#ifdef PF
      case 'm': num_particles = atoi(optarg); break;
      case 'W': reverseTrees = ! reverseTrees; break;
#endif
      default: break;
    }
  }

  return 0;
}

// ===========================================================================

#ifdef USE_MPI

double lik_master(Forest& trees, const vector<double>& vars,
                  int SImodel, int vflag, int rescale) 
{
  double fx(0.0);
  double f(0.0);
  int id(0);
  MPI::Status status;
  int flag;
  int src;
  int idle(0);
  int k(0);
  size_t j(0);
  for (k = 1; k < mpi_ntasks && j < trees.size(); ++k) {
    MPI::COMM_WORLD.Send(&j,1,MPI::INT,k,1);
    MPI::COMM_WORLD.Send(vars.data(),vars.size(),MPI::DOUBLE,k,1);
    // if (vflag) fprintf(stderr,"Sending tree %d to process %d.\n",(int) j,k);
    ++j;
  }
  while (1) {
    flag = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE,MPI::ANY_TAG,status);
    if (flag) {
      src = status.Get_source();
      MPI::COMM_WORLD.Recv(&id,1,MPI::INT,MPI::ANY_SOURCE,MPI::ANY_TAG,status);
      MPI::COMM_WORLD.Recv(&f,1,MPI::DOUBLE,src,MPI::ANY_TAG,status);
      // if (vflag) fprintf(stderr,"Receiving tree %d from process %d.\n",id,src);
      if (j < trees.size()) {
        MPI::COMM_WORLD.Send(&j,1,MPI::INT,src,1);
        MPI::COMM_WORLD.Send(vars.data(),vars.size(),MPI::DOUBLE,src,1);
        ++j;
      } else {
        ++idle;
      }
      fx += f; 
      if (idle == mpi_ntasks-1) break;
    }
  }
  return fx;
}

// ===========================================================================

void lik_slave(Forest& trees, vector<double>& vars, 
                  int SImodel, int vflag, int rescale, int survival) {
  double f(0.0);
  int id(0);
  MPI::Status status;
  while (1) {
    MPI::COMM_WORLD.Recv(&id,1,MPI::INT,0,MPI::ANY_TAG,status);
    if (status.Get_tag() == 0) break;
    else {
      MPI::COMM_WORLD.Recv(vars.data(),vars.size(),MPI::DOUBLE,0,MPI::ANY_TAG,status);
      last_time = clock();
#ifdef PF
      f = -INFINITY;
//       f = pfLikelihood(vars,vars,num_particles,rng);
#else
      f = trees.at(id)->likelihood(vars,SImodel,vflag,rescale,survival);
#endif
      calc_time += clock()-last_time;
      MPI::COMM_WORLD.Send(&id,1,MPI::INT,0,2);
      MPI::COMM_WORLD.Send(&f,1,MPI::DOUBLE,0,2);
    }
  }
  fprintf(stderr,"Slave %d done.\n",mpi_rank);
}


#endif
