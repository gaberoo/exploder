#include "TreeSwarm.h"

double likelihood(const PSO::Point& x, void* params) {
  PSOParams par = *(PSOParams*) params;

  double fx;
  int    N    = (int) round(x[4]); 
  double beta = x[0];
  double mu   = (par.lockedRatio >= 0.0) ? x[2]*(1./par.lockedRatio-1.) : x[1];
  double psi  = x[2];
  double rho  = x[3];
  int vflag = (par.vflag > 0) ? par.vflag-1 : vflag;

  vector<double> vars(x);
  vars[1] = mu;

  if ((N <= par.maxExtant && par.SImodel) || beta <= 0.0 || mu < 0.0 || psi < 0.0 
      || rho < 0.0 || rho > 1.0) {
    fx = -INFINITY;
  } else {
    if (! par.SImodel) vars[4] = 0.0;
    if (par.tree >= 0)
      fx = par.trees->likelihood(par.tree,vars,par.SImodel,vflag,par.rescale);
    else 
      fx = par.trees->jointLikelihood(vars,par.SImodel,vflag,par.rescale);
  }

  return fx;
}


#ifdef USE_MPI
/*
void PSO::TreeSwarm::evaluate_master(int vflag) {
  int id(0);
  int j(0);
  int tree(0);
  MPI::Status status;
  int flag;
  int src;
  int idle(0);
  PSOParams* par = (PSOParams*) p->evalParams;
  const int numTrees(par->trees->size());

  queue< pair<int,int> > evalQueue;  // pair<Particle,Tree>
  vector<int> treesWaiting(swarm.size(),numTrees);
  for (int i(0); i < swarmSize; ++i) {
    for (int t(0); t < numTrees; ++t) {
      evalQueue.push(make_pair(i,t));
    }
  }

  double f(-INFINITY);
  double mf;
  vector< vector<double> > ff(swarm.size(),vector<double>(par->trees->size()));
  int ff_cnt;

  if (vflag) cerr << "Sending particles to slaves..." << endl;
  // initialize slaves
  for (int k(1); k < mpi_ntasks && ! evalQueue.empty(); ++k) {
    j = evalQueue.front().first;
    tree = evalQueue.front().second;
    evalQueue.pop();
    if (vflag) fprintf(stderr,"Sending particle %d/%d to process %d.\n",j,tree,k);
    MPI::COMM_WORLD.Send(&j,1,MPI::INT,k,1);
    MPI::COMM_WORLD.Send(swarm[j]->pos(),numParams,MPI::DOUBLE,k,1);
    MPI::COMM_WORLD.Send(&tree,1,MPI::INT,k,1);
  }

  while (1) {
    flag = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE,MPI::ANY_TAG,status);
    if (flag) {
      // get function value
      src = status.Get_source();
      MPI::COMM_WORLD.Recv(&id,1,MPI::INT,src,MPI::ANY_TAG,status);
      MPI::COMM_WORLD.Recv(&tree,1,MPI::INT,src,MPI::ANY_TAG,status);
      MPI::COMM_WORLD.Recv(&f,1,MPI::DOUBLE,src,MPI::ANY_TAG,status);
      if (vflag) fprintf(stderr,"Receiving particle %d/%d from process %d.\n",id,tree,src);

      ff[id][tree] = f;

      --treesWaiting[id];
      if (treesWaiting[id] == 0) {
        cerr << "calc'ing..." << endl;
        f = 0.0;
        mf = 0.0;
        ff_cnt = 0;
        for (int t(0); t < numTrees; ++t) {
          if (isfinite(ff[id][t])) { f += ff[id][t]; ++ff_cnt; }
        }
        f /= ff_cnt;
        for (int t(0); t < numTrees; ++t) mf += exp(ff[id][t]-f);
        f = log(mf)+f-log(ff_cnt);

        // update particle information
        swarm[id]->value = f;
        if (f >= swarm[id]->bestValue) {
          swarm[id]->bestPosition = swarm[id]->position;
          swarm[id]->bestValue = f;
        }
        ++numEvals;

        // check for new best value
        if (f >= bestVal) {
          bestPos = swarm[id]->position;
          bestVal = f;
          bestParticle = id;
        }
      }

      // send new work to slave
      if (! evalQueue.empty()) {
        j = evalQueue.front().first;
        tree = evalQueue.front().second;
        evalQueue.pop();
        if (vflag) fprintf(stderr,"Sending particle %d/%d to process %d.\n",j,tree,src);
        MPI::COMM_WORLD.Send(&j,1,MPI::INT,src,1);
        MPI::COMM_WORLD.Send(swarm[j]->position.data(),numParams,MPI::DOUBLE,src,1);
        MPI::COMM_WORLD.Send(&tree,1,MPI::INT,src,1);
      } else {
        if (vflag) fprintf(stderr,"Sending done signal to process %d.\n",src);
        MPI::COMM_WORLD.Send(0,0,MPI::INT,src,0);
        ++idle;
      }

      if (idle == mpi_ntasks-1) break;
    }
  }
}
*/
// ===========================================================================

void PSO::TreeSwarm::run_master(int numIt, int vflag, ostream* out, ostream* hist) {
  int id(0);
  int j(0);
  int tree(0);
  MPI::Status status;
  int flag;
  int src;
  int idle(0);
  int iter(0);
  PSOParams* par = (PSOParams*) p->evalParams;
  const int numTrees(par->trees->size());

  queue< pair<int,int> > evalQueue;  // pair<Particle,Tree>
  vector<int> treesWaiting(swarm.size(),numTrees);
  for (int i(0); i < swarmSize; ++i) {
    if (numIt > 0) {
      updateVelocity(i);
      updatePosition(i);
    }
    for (int t(0); t < numTrees; ++t) {
      evalQueue.push(make_pair(i,t));
    }
  }

  double f(-INFINITY);
  double mf;
  vector< vector<double> > ff(swarm.size(),vector<double>(par->trees->size()));
  int ff_cnt;

  if (vflag) cerr << "Sending particles to slaves..." << endl;
  // initialize slaves
  for (int k(1); k < mpi_ntasks && (iter < numIt || numIt < 0); ++k) {
    j = evalQueue.front().first;
    tree = evalQueue.front().second;
    evalQueue.pop();
    if (vflag) cerr << j << " " << (*swarm[j]) << endl;
    if (vflag) fprintf(stderr,"Sending particle %d/%d to process %d.\n",j,tree,k);
    MPI::COMM_WORLD.Send(&j,1,MPI::INT,k,1);
    MPI::COMM_WORLD.Send(swarm[j]->position.data(),numParams,MPI::DOUBLE,k,1);
    MPI::COMM_WORLD.Send(&tree,1,MPI::INT,k,1);
    ++iter;
  }

  while (1) {
    flag = MPI::COMM_WORLD.Iprobe(MPI::ANY_SOURCE,MPI::ANY_TAG,status);
    if (flag) {
      // get function value
      src = status.Get_source();
      MPI::COMM_WORLD.Recv(&id,1,MPI::INT,src,MPI::ANY_TAG,status);
      if (vflag) fprintf(stderr,"Receiving particle %d from process %d.\n",id,src);
      MPI::COMM_WORLD.Recv(&tree,1,MPI::INT,src,MPI::ANY_TAG,status);
      MPI::COMM_WORLD.Recv(&f,1,MPI::DOUBLE,src,MPI::ANY_TAG,status);

      ff[id][tree] = f;

      if (hist != NULL) {
        *hist << id << " " << tree << " "  << f << " "
              << swarm[id]->position << endl;
      }

      --treesWaiting[id];
      if (treesWaiting[id] == 0) {
        /*
        f = 0.0;
        mf = 0.0;
        ff_cnt = 0;
        for (int t(0); t < numTrees; ++t) {
          if (isfinite(ff[id][t])) { f += ff[id][t]; ++ff_cnt; }
        }
        f /= ff_cnt;
        for (int t(0); t < numTrees; ++t) {
          if (isfinite(ff[id][t])) mf += exp(ff[id][t]-f);
        }
        if (ff_cnt > 0) f = log(mf)+f-log(ff_cnt);
        else f = -INFINITY;
        */
        vector<double>::const_iterator fi(ff[id].begin());
        f = 0.0;
        do f += *fi++; while (fi != ff[id].end());
        // f = accumulate(ff[id].begin(),ff[id].end(),0.0);
        f /= numTrees;

        // update particle information
        swarm[id]->value = f;
        if (f >= swarm[id]->bestValue) {
          swarm[id]->bestPosition = swarm[id]->position;
          swarm[id]->bestValue = f;
        }
        ++numEvals;

        // check for new best value
        if (f >= bestVal) {
          bestPos = swarm[id]->position;
          bestVal = f;
          bestParticle = id;
          if (out != NULL) {
            *out << numEvals << " " << bestVal << " ";
            for (size_t j(0); j < bestPos.size(); ++j) *out << bestPos[j] << " ";
            *out << endl;
          }
        }

        if (numIt > 0) {
          // update velocity and position
          updateVelocity(id);
          updatePosition(id);
          for (int t(0); t < numTrees; ++t) evalQueue.push(make_pair(id,t));
          treesWaiting[id] += numTrees;
        }
      }

      // send new work to slave
      if ((iter < numIt || numIt < 0) && ! evalQueue.empty()) {
        j = evalQueue.front().first;
        tree = evalQueue.front().second;
        evalQueue.pop();
        if (vflag) fprintf(stderr,"Sending particle %d/%d to process %d.\n",j,tree,src);
        MPI::COMM_WORLD.Send(&j,1,MPI::INT,src,1);
        MPI::COMM_WORLD.Send(swarm[j]->position.data(),numParams,MPI::DOUBLE,src,1);
        MPI::COMM_WORLD.Send(&tree,1,MPI::INT,src,1);
        ++iter;
      } else {
        ++idle;
        if (vflag) fprintf(stderr,"Sending done signal to process %d. %d running processes left.\n",src,mpi_ntasks-idle);
        MPI::COMM_WORLD.Send(0,0,MPI::INT,src,0);
      }

      if (idle == mpi_ntasks-1) break;
    }
  }
}

// ===========================================================================

void PSO::TreeSwarm::evaluate_slave() {
  double f(-INFINITY);
  int id(0);
  int tree(-1);
  Point position(numParams);
  MPI::Status status;
  PSOParams* par = (PSOParams*) p->evalParams;
  while (1) {
    MPI::COMM_WORLD.Recv(&id,1,MPI::INT,0,MPI::ANY_TAG,status);
    if (status.Get_tag() == 0) break;
    MPI::COMM_WORLD.Recv(position.data(),numParams,MPI::DOUBLE,0,MPI::ANY_TAG,status);
    MPI::COMM_WORLD.Recv(&tree,1,MPI::INT,0,MPI::ANY_TAG,status);
    par->tree = tree;
    f = p->evalFunc(position,p->evalParams);
    MPI::COMM_WORLD.Send(&id,1,MPI::INT,0,2);
    MPI::COMM_WORLD.Send(&tree,1,MPI::INT,0,2);
    MPI::COMM_WORLD.Send(&f,1,MPI::DOUBLE,0,2);
  }
}
#endif
