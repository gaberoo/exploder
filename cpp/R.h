#ifndef __FAKE_R__
#define __FAKE_R__

#define Rprintf(format, args...) fprintf(stderr, format , ## args)
#define R_PosInf INFINITY
#define R_NegInf -R_PosInf

#endif // __FAKE_R__
