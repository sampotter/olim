#ifndef __UPDATE_RULES_LINE_UPDATES_HPP__
#define __UPDATE_RULES_LINE_UPDATES_HPP__

#include "cost_funcs.hpp"

namespace updates {

template <cost_func F, int d>
struct line_bv {
  inline double operator()(double u0, double s, double s0, double h) const {
    if (F == RHR) {
      (void) s0;
#if PRINT_UPDATES
      double u = u0 + s*h*sqrt(d);
      printf("line<1>(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",u0,s,s0,h,u);
      return u;
#else
      return u0 + s*h*sqrt(d);
#endif
    } else {
      (void) s0;
#if PRINT_UPDATES
      double u = u0 + (s + s0)*h*sqrt(d)/2;
      printf("line<1>(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",u0,s,s0,h,u);
      return u;
#else
      return u0 + (s + s0)*h*sqrt(d)/2;
#endif
    }
  }
};

template <cost_func F, int n>
struct line {
  double operator()(double const * p0, double u0, double s, double s0,
                    double h) const {
    (void) s0;
    if (F == RHR) {
#if PRINT_UPDATES
      double u = u0 + s*h*norm2<n>(p0);
      printf("line<%d>(p0 = {", n);
      for (int i = 0; i < n - 1; ++i) printf("%g, ", p0[i]);
      printf("%g}, ", p0[n - 1]);
      printf("u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",u0,s,s0,h,u);
      return u;
#else
      return u0 + s*h*norm2<n>(p0);
#endif
    } else {
#if PRINT_UPDATES
      double u = u0 + (s + s0)*h*norm2<n>(p0)/2;
      printf("line<%d>(p0 = {", n);
      for (int i = 0; i < n - 1; ++i) printf("%g, ", p0[i]);
      printf("%g}, ", p0[n - 1]);
      printf("u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",u0,s,s0,h,u);
      return u;
#else
      return u0 + (s + s0)*h*norm2<n>(p0)/2;
#endif
    }
  }
};

}

#endif // __UPDATE_RULES_LINE_UPDATES_HPP__
