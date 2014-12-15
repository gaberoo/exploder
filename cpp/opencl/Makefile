include ../Make.inc

# Yosemite broke the OpenCL framework for MacPorts GCC
CC = clang
LDFLAGS += -L/opt/local/lib/gcc49 

CFLAGS += -I. -I.. $(CL_CFLAGS)
LDFLAGS += $(CL_LDFLAGS)

SRCS = $(wildcard *.c)
OBJS = $(SRCS:.c=.o)

TEST_SRC = tests/$(wildcard *.c)
TEST_OBJ = $(TEST_SRC:.c=.o)
TESTS = $(TEST_SRC:.c=)

EXPOMV = $(shell ls ../shared/*.{c,cpp,f})
EXPOMV_LIB = ../shared/expomv.a

clean:
	rm -rf $(OBJS)

distclean: clean
	rm -f expocl.a
	rm -f $(TEST_OBJ)

tests: $(TEST_OBJ)
	make $(TESTS)

lib: expocl.a

expocl.a: $(OBJS)
	ar rcs expocl.a $(OBJS)

../shared/expomv.a: $(EXPOMV)
	cd ..; make expomv.a; cd opencl;

../Tree.o: ../Tree.h ../Tree.cpp
	cd ..; make Tree.o; cd opencl;

../expoTree.o: ../expoTree.h ../expoTree.cpp
	cd ..; make expoTree.o; cd opencl;

##############################################################################

tests/test_compile: tests/test_compile.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

test_expocl: test_expocl.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o test_expocl test_expocl.c expocl.a $(EXPOMV_LIB) \
		$(LDFLAGS) -lgsl 

test_expocl_taylor: test_expocl_taylor.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

tests/test_expocl_inf_norm: tests/test_expocl_inf_norm.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

tests/test_expocl_matvec: tests/test_expocl_matvec.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

tests/test_expocl_normam: tests/test_expocl_normam.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

tests/test_expocl_expmv: tests/test_expocl_expmv.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

tests/test_expocl_norms: tests/test_expocl_norms.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

tests/test_expocl_expotree: tests/test_expocl_expotree.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

test_expocl_stdeg: test_expocl_stdeg.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

test_expocl_shifts: test_expocl_shifts.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

test_expocl_neg_vals: test_expocl_neg_vals.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

tests/test_expocl_int_ival: tests/test_expocl_int_ival.c expocl.a $(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ expocl.a $(EXPOMV_LIB) $(LDFLAGS) -lgsl 

test_cl_inf_norm: test_cl_inf_norm.c util.c
	$(CC) $(CFLAGS) -o $@ $^

test_cl_tridiag: test_cl_tridiag.c util.c
	$(CC) $(CFLAGS) -o $@ $^ -lgsl

tests/test_expocl_tridiag_index: tests/test_expocl_tridiag_index.c expocl.a \
	$(EXPOMV_LIB)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lgsl 

##############################################################################

.c.o: $<
	$(CC) $(CFLAGS) -c $<

.f.o: $<
	$(FC) $(FFLAGS) -c $<

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

