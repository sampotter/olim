#ifndef __FMM_MEX_HPP__
#define __FMM_MEX_HPP__

#include <cstddef>

#include "typedefs.hpp"

enum class marcher_type {
    basic
};

extern "C"
void fmm_mex(double * out, bool * in, size_t M, size_t N, double h,
			 speed_func F, marcher_type type);

#endif // __FMM_MEX_HPP__
