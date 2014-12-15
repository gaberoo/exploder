include ../Make.inc

LDFLAGS += -lgsl
CFLAGS += -I.. -I../shared

BLAS = -framework Accelerate

##############################################################################

all: distclean lib

clean:
	@echo "Cleaning build files..."
	@rm -rf *.o;

distclean: clean

##############################################################################

lib:
	cd ..; make shared/expomv.a; cd tests;

##############################################################################

test_find_ival: test_find_ival.c ../shared/find_intervals.c
	$(CC) $(CFLAGS) -o $@ $^

##############################################################################

makeMat: makeMat.c lib
	$(CC) $(CFLAGS) -o makeMat makeMat.c ../expomv.a -framework Accelerate

test_one_norm: test_one_norm.c lib
	$(CC) $(CFLAGS) -o $@ $< ../expomv.a -framework Accelerate $(LDFLAGS)

testNorm: testNorm.c lib
	$(CC) $(CFLAGS) -o testNorm testNorm.c ../expomv.a $(LDFLAGS)

testExpmv: testExpmv.c lib dgchbv.f
	$(FC) $(FFLAGS) -c dgchbv.f
	$(CC) $(CFLAGS) -o testExpmv testExpmv.c ../expomv.a dgchbv.o $(LDFLAGS)

testCopy: testCopy.c lib
	$(CC) $(CFLAGS) -o testCopy testCopy.c ../expomv.a $(LDFLAGS)

testDlange: testDlange.c
	$(CC) $(CFLAGS) -o testDlange testDlange.c $(LDFLAGS)

full_prob: full_prob.cpp ../expoTree.o ../expomv.a ../Tree.o
	$(CPP) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lgsl

test_json: test_json.cpp
	cd ..; make json.o; cd tests;
	$(CPP) $(CFLAGS) -o test_json test_json.cpp ../json.o $(LDFLAGS)

test_bdsky: test_bdsky.c lib
	$(CC) $(CFLAGS) -o $@ $< ../shared/expomv.a

# EXT ########################################################################

makeMat_ext: makeMat_ext.c lib
	$(CC) $(CFLAGS) -o makeMat_ext makeMat_ext.c ../expomv.a -framework Accelerate

test_ext: test_ext.cpp lib
	cd ..; make expoTree.o; make expoTreeExt.o; cd tests;
	$(CPP) $(CFLAGS) -o test_ext test_ext.cpp ../expomv.a \
		../expoTree.o ../expoTreeExt.o \
		-lgsl -framework Accelerate

##############################################################################

test_sir_index: test_sir_index.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

test_sir_mat: test_sir_mat.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

test_sir_init_mat: test_sir_init_mat.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

test_sir_matfun: test_sir_matfun.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

test_sir_one_norm: test_sir_one_norm.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

test_sir_inf_norm: test_sir_inf_norm.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

test_sir_trans: test_sir_trans.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

test_sir_expmv: test_sir_expmv.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

test_sir_print_mat: test_sir_print_mat.c lib
	$(CC) $(CFLAGS) $< -o $@ ../expomv.a $(BLAS)

