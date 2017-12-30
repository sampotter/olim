#ifndef __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.defs.hpp"
#include "common.macros.hpp"
#include "olim_util.hpp"
#include "update_rules.tetra_updates.util.hpp"

template <class speed_est, char degree>
template <char p0, char p1, char p2>
double update_rules::tetra_updates<speed_est, degree>::tetra(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  ffvec<p0>, ffvec<p1>, ffvec<p2>) const
{
  return tetra_impl<p0, p1, p2>(
    u0, u1, u2, s, s0, s1, s2, h,
    ffvec<p0> {}, ffvec<p1> {}, ffvec<p2> {},
    std::integral_constant<char, degree> {});
}

template <class speed_est, char degree>
template <char p0, char p1, char p2>
double update_rules::tetra_updates<speed_est, degree>::tetra_impl(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  ffvec<p0>, ffvec<p1>, ffvec<p2>,
  std::integral_constant<char, 0>) const
{
  (void) u0;
  (void) u1;
  (void) u2;
  (void) s;
  (void) s0;
  (void) s1;
  (void) s2;
  (void) h;
  return INF(double);
}

template <class speed_est, char degree>
template <char p0, char p1, char p2>
double update_rules::tetra_updates<speed_est, degree>::tetra_impl(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  ffvec<p0>, ffvec<p1>, ffvec<p2>,
  std::integral_constant<char, 1>) const
{
  (void) u0;
  (void) u1;
  (void) u2;
  (void) s;
  (void) s0;
  (void) s1;
  (void) s2;
  (void) h;
  return INF(double);
}

#endif // __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__
