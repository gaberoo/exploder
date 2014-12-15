#include "expmv.h"
#include "expo_type.h"

void taylor_length(double tt, const double* tm, const void* pars,
                     double* s, int* tcol) 
{
  const expo_type* expo = (expo_type*) pars;

  /* calculate costs */
  double cost = R_PosInf;
  double ccost = 0.0;

  double Ck = 0.0;
  *tcol = 0;
  *s = 0.0;

  int i, j;
  for (i = 0; i < expo->m_max; ++i) {
    ccost = R_PosInf;
    for (j = 0; j < expo->p_max-1; ++j) {
      Ck = ceil(fabs(tt)*tm[j*(expo->m_max)+i])*(i+1.);
      if (Ck == 0.0) Ck = R_PosInf;
      if (Ck < ccost) ccost = Ck;
    }
    if (ccost < cost) {
      cost = ccost;
      *tcol = i+1;
    }
  }

  /* get length of Taylor explansion */
  *s = fmax(cost/(*tcol),1.0);
}

