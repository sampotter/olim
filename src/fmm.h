#ifndef __FMM_H__
#define __FMM_H__

#include "typedefs.h"

extern "C"
void fmm(double * out, bool * in, int M, int N, double h, double * pvalues,
         marcher_type type);

extern "C"
void fmm3d(double * out, bool * in, int* dims, double h, double * S,
           marcher_type type);

#endif // __FMM_H__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
