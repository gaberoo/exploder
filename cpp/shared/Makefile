include ../Make.inc

CFLAGS += -I. -I..

C_SRCS = $(wildcard *.c)
F_SRCS = $(wildcard *.f)

OBJS  = $(C_SRCS:.c=.o)
OBJS += $(F_SRCS:.f=.o)

clean:
	rm -rf $(OBJS)

distclean: clean
	rm -f expomv.a

lib: expomv.a

expomv.a: $(OBJS)
	ar rcs expomv.a $(OBJS)

##############################################################################

.c.o: $<
	$(CC) $(CFLAGS) -c $<

.f.o: $<
	$(FC) $(FFLAGS) -c $<

.cpp.o: $<
	$(CPP) $(CPPFLAGS) -c $<

