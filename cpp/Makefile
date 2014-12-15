include Make.inc

ifdef DEBUG
CFLAGS += -DDEBUG
endif

CFLAGS += -isystem ./include

##############################################################################

PROGS = testLik calc_likelihood sim_trees dream

##############################################################################

C_SRCS = $(wildcard *.c)
F_SRCS = $(wildcard *.c)
CPP_SRCS = $(wildcard *.cpp)

OBJS = $(C_SRCS:.c=.o)
OBJS += $(F_SRCS:.f=.o)
OBJS += $(CPP_SRCS:.cpp=.o)

DLACN1 = dlacn1.a
EXPOMV = shared/expomv.a
EXPOCL = opencl/expocl.a
MCTRAJ = mctraj/mctraj.a

PSO_DIR = pso
PSO = $(PSO_DIR)/libpso.a
MPIPSO = $(PSO_DIR)/libpso_mpi.a

SUBPLEX_DIR = subplex
SUBPLEX = $(SUBPLEX_DIR)/subplex.a

##############################################################################

all:

clean:
	rm -f *.o;
	rm -f obj/*.o; 
	rm -f $(PROGS);

distclean:
	make clean;
	cd pso; make distclean;
	cd mctraj; make distclean;
	cd opencl; make distclean;
	cd shared; make distclean;

version:
	$(CC) --version

##############################################################################

$(EXPOMV): shared/*.h shared/*.c shared/*.f shared/*.cpp
	cd shared; make expomv.a

$(EXPOCL): opencl/*.c
	cd opencl; make expocl.a

$(PSO): $(PSO_DIR)/*.h $(PSO_DIR)/*.cpp
	cd pso; make; cd ..;

$(MCTRAJ): mctraj/*.h mctraj/*.cpp
	cd mctraj; make mctraj.a

$(SUBPLEX):
	cd $(SUBPLEX_DIR); make subplex.a; cd ..;

r_sim_trees.o: sim_trees.cpp
	$(CPP) $(CPPFLAGS) -DRLIB -c -o $@ $^

##############################################################################

testLik: testLik.o expoTree.o expomv.o $(EXPOMV)
	$(LD) -o $@ $^ -lgsl $(LDFLAGS)

sim_trees: sim_trees.cpp $(EXPOMV)
	$(CPP) -o $@ $^ -lmygsl -lgsl $(LDFLAGS)

sim_epi: sim_epi.cpp
	$(CPP) $(CPPFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

calc_likelihood: calc_likelihood.cpp expoTree.o \
	expoTreeSIR.o $(EXPOMV) Tree.o
	$(LD) -o $@ $^ $(LDFLAGS) -lgsl

full_prob: full_prob.cpp expoTree.o \
	expoTreeSIR.o expoTreePSave.o $(EXPOMV) Tree.o
	$(LD) -o $@ $^ $(LDFLAGS) -lgsl

calc_likelihood_cl: calc_likelihood.cpp expoTree.cpp \
	expoTreeSIR.cpp Tree.cpp $(EXPOMV) $(EXPOCL)
	$(CPP) -DOPENCL -o $@ $^ $(LDFLAGS) -lgsl $(CL_LDFLAGS)

optim: optim.cpp expoTree.o expoTreeSIR.o $(EXPOMV) Tree.o $(SUBPLEX)
	$(LD) -o $@ $^ $(LDFLAGS) -lgsl

optim_cl: optim.cpp expoTree.cpp \
	expoTreeSIR.cpp Tree.cpp $(EXPOMV) $(EXPOCL) $(SUBPLEX)
	$(CPP) -DOPENCL -o $@ $^ $(LDFLAGS) -lgsl $(CL_LDFLAGS)

optim_pso: optim_pso.o expomv.o expoTree.o $(EXPOMV) TreeSwarm.o $(PSO)
	$(CPP) $(CPPFLAGS) -I$(PSO_DIR) -o $@ $^ -lgsl $(LDFLAGS)

optim_pso_inf: optim_pso.cpp expomv.o expoTree.o $(EXPOMV) $(PSO) TreeSwarm.cpp
	$(CPP) $(CPPFLAGS) -I$(PSO_DIR) -DINF_N -o $@ $^ -lgsl $(LDFLAGS)

optim_pso_mpi: optim_pso.cpp expomv.o expoTree.o $(EXPOMV) TreeSwarm.cpp $(MPIPSO)
	$(MPICPP) $(CPPFLAGS) -DUSE_MPI -I$(PSO_DIR) -o $@ $^ -lgsl $(MPILDFLAGS)

mcmc: mcmc.cpp expomv.o expoTree.o $(EXPOMV)
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

mcmc_mpi: mcmc.cpp expomv.o expoTree.o $(EXPOMV)
	$(MPICPP) -DUSE_MPI $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS) $(MPILDFLAGS)

mcmc_inf: expomv.o expoTree.o $(EXPOMV) mcmc_inf.o
	$(LD) -o $@ $^ -lgsl $(LDFLAGS)

dream: dream.cpp gelman_rubin.o expoTree.o expoTreeSIR.o $(EXPOMV) \
	Forest.h Tree.h Tree.cpp
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

dream-shifts: dream-shifts.cpp gelman_rubin.o expoTree.o expoTreeSIR.o $(EXPOMV) \
	Forest.h Tree.h Tree.cpp
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

dream_mpi: dream.cpp gelman_rubin.o expoTree.o expoTreeSIR.o $(EXPOMV) \
	Forest.h Tree.h Tree.cpp
	$(MPICPP) -DUSE_MPI $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS) $(MPILDFLAGS)

dream_pf: dream.cpp gelman_rubin.o expoTree.o expoTreeSIR.o $(EXPOMV) \
	Forest.h Tree.h Tree.cpp $(MCTRAJ)
	$(CPP) -DPF -fopenmp $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

dream_pf_mpi: dream.cpp gelman_rubin.o expoTree.o expoTreeSIR.o $(EXPOMV) \
	Forest.h Tree.h Tree.cpp $(MCTRAJ)
	$(MPICPP) -DUSE_MPI -DPF $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS) $(MPILDFLAGS)

testGelman: testGelman.cpp gelman_rubin.o
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

det_epi: det_epi.cpp
	$(CPP) $(CFLAGS) -o $@ $^ -lgsl $(LDFLAGS)

##############################################################################

calc_likelihood_ext: calc_likelihood_ext.cpp expoTreeExt.o expoTree.o \
	$(EXPOMV) Tree.o expoTreeSIR.o
	$(CPP) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lgsl

##############################################################################

.c.o:; $(CC) $(CFLAGS) -c $<

.f.o:; $(FC) $(FFLAGS) -c $<

.cpp.o:; $(CPP) $(CPPFLAGS) -c $<


