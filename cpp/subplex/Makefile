#Makefile for building subplex

FORTRAN = gfortran
OPTS    = -O3
LOADER  = gfortran

OBJECTS = demo.o

SUBPLEXLIB = subplex.a

#Uncomment next two lines to link with system BLAS
SUBBLASLIB =
SYSBLASLIB = -framework Accelerate
#Uncomment next two lines to link with subplex BLAS provided with distribution
#SUBBLASLIB = subblas.a
#SYSBLASLIB =

SUBPLEXOBJECTS = calcc.o dist.o evalf.o fstats.o newpt.o order.o \
                 partx.o setstp.o simplx.o sortd.o start.o subopt.o \
                 subplx.o
SUBBLASOBJECTS = dasum.o daxpy.o dcopy.o dscal.o

demo: $(OBJECTS) $(SUBPLEXLIB) $(SUBBLASLIB)
	$(LOADER) $(OBJECTS) $(SUBPLEXLIB) $(SUBBLASLIB) $(SYSBLASLIB) -o demo

subplex.a: $(SUBPLEXOBJECTS)
	ar cr subplex.a $(SUBPLEXOBJECTS)
	ranlib subplex.a

subblas.a: $(SUBBLASOBJECTS)
	ar cr subblas.a $(SUBBLASOBJECTS)
	ranlib subblas.a

.f.o: ; $(FORTRAN) $(OPTS) -c $<
