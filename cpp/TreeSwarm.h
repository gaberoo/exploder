#ifndef __TREE_SWARM__
#define __TREE_SWARM__

#include "pso/Point.h"
#include "pso/Swarm.h"

#include "expoTree.h"
#include "Forest.h"

#ifdef USE_MPI
#include "mpi.h"
#endif

// ===========================================================================

typedef struct {
  int extant;
  int maxExtant;
  int SImodel;
  int vflag;
  int rescale;
  double lockedRatio;
  Forest* trees;
  int tree;
} PSOParams;

double likelihood(const PSO::Point& x, void* params);

// ===========================================================================

namespace PSO {
  class TreeSwarm : public Swarm {
    public:
      TreeSwarm() : Swarm() { p->evalFunc = &likelihood; }
      TreeSwarm(int size, int np, Parameters* pars) : Swarm(size,np,pars) { p->evalFunc = &likelihood; }
      TreeSwarm(int size, int np, Parameters* pars, PSOParams* evalPars) 
        : Swarm(size,np,pars)
      { 
        p->evalFunc = &likelihood; 
        p->evalParams = evalPars; 
      }
      virtual ~TreeSwarm() {}

      void setEvalParams(PSOParams* pars) { p->evalParams = pars; }

#ifdef USE_MPI
      void evaluate_slave();
      void run_master(int numInt, int vflag, ostream* out, ostream* hist);
#endif
  };
};

#endif // __TREE_SWARM__
