#ifndef __UPDATE_RULES_LINE_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_LINE_UPDATES_IMPL_HPP__

#include "common.defs.hpp"
#include "update_rules.utils.hpp"

#include <src/config.hpp>

#if PRINT_UPDATES
#    include <cstdio>
#endif

namespace update_rules {

template <>
inline double
rhr_line_updates::line<1>(double u0, double s, double s0, double h) const {
  (void) s0;
#if PRINT_UPDATES
  double tmp = u0 + s*h;
  printf("line<1>(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",
         u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + s*h;
#endif
}

template <>
inline double
rhr_line_updates::line<2>(double u0, double s, double s0, double h) const {
  (void) s0;
#if PRINT_UPDATES
  double tmp = u0 + s*h*sqrt2;
  printf("line<2>(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",
         u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + s*h*sqrt2;
#endif
}

template <>
inline double
rhr_line_updates::line<3>(double u0, double s, double s0, double h) const {
  (void) s0;
#if PRINT_UPDATES
  double tmp = u0 + s*h*sqrt3;
  printf("line<3>(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",
         u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + s*h*sqrt3;
#endif
}

template <>
inline double
rhr_line_updates::line<2>(double const * p0, double u0, double s, double s0,
                          double h) const {
  (void) s0;
#if PRINT_UPDATES
#  error Not implemented
#else
  return u0 + s*h*norm2<2>(p0);
#endif
}

template <>
inline double
rhr_line_updates::line<3>(double const * p0, double u0, double s, double s0,
                          double h) const {
  (void) s0;
#if PRINT_UPDATES
#  error Not implemented
#else
  return u0 + s*h*norm2<3>(p0);
#endif
}

template <>
inline double
mp_line_updates::line<1>(double u0, double s, double s0, double h) const {
#if PRINT_UPDATES
  double tmp = u0 + (s + s0)*h/2;
  printf("line<1>(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",
         u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + (s + s0)*h/2;
#endif
}

template <>
inline double
mp_line_updates::line<2>(double u0, double s, double s0, double h) const {
#if PRINT_UPDATES
  double tmp = u0 + (s + s0)*h*sqrt2/2;
  printf("line<2>(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n",
         u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + (s + s0)*h*sqrt2/2;
#endif
}

template <>
inline double
mp_line_updates::line<3>(double u0, double s, double s0, double h) const {
#if PRINT_UPDATES
  double tmp = u0 + (s + s0)*h*sqrt3/2;
  printf("line<3>(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n", u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + (s + s0)*h*sqrt3/2;
#endif
}

template <>
inline double
mp_line_updates::line<2>(double const * p0, double u0, double s, double s0,
                         double h) const {
#if PRINT_UPDATES
#  error Not implemented
#else
  return u0 + (s + s0)*h*norm2<2>(p0)/2;
#endif
}

template <>
inline double
mp_line_updates::line<3>(double const * p0, double u0, double s, double s0,
                         double h) const {
#if PRINT_UPDATES
#  error Not implemented
#else
  return u0 + (s + s0)*h*norm2<3>(p0)/2;
#endif
}

}

#endif // __UPDATE_RULES_LINE_UPDATES_IMPL_HPP__
