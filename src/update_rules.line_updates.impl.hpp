#ifndef __UPDATE_RULES_LINE_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_LINE_UPDATES_IMPL_HPP__

#include "update_rules.line_updates.hpp"

#include "common.defs.hpp"

#include <src/config.hpp>

#if PRINT_UPDATES
#    include <cstdio>
#endif

template <>
inline double
update_rules::line_updates::line<1>(double u0, double s, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + s*h;
  printf("line1(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h;
#endif
}

template <>
inline double
update_rules::line_updates::line<1>(double u0, double s, double s0, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + (s + s0)*h/2;
  printf("line1(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n", u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + (s + s0)*h/2;
#endif
}

template <>
inline double
update_rules::line_updates::line<2>(double u0, double s, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + s*h*sqrt2;
  printf("line2(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h*sqrt2;
#endif
}

template <>
inline double
update_rules::line_updates::line<2>(double u0, double s, double s0, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + (s + s0)*h*sqrt2/2;
  printf("line2(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n", u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + (s + s0)*h*sqrt2/2;
#endif
}

template <>
inline double
update_rules::line_updates::line<3>(double u0, double s, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + s*h*sqrt3;
  printf("line3(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h*sqrt3;
#endif
}

template <>
inline double
update_rules::line_updates::line<3>(double u0, double s, double s0, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + (s + s0)*h*sqrt3/2;
  printf("line3(u0 = %g, s = %g, s0 = %g, h = %g) -> %g\n", u0, s, s0, h, tmp);
  return tmp;
#else
  return u0 + (s + s0)*h*sqrt3/2;
#endif
}

#endif // __UPDATE_RULES_LINE_UPDATES_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
