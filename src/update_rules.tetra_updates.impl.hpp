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
#include "numopt.hpp"
#include "olim_util.hpp"
#include "update_rules.tetra_updates.util.hpp"

template <class derived>
template <char p0, char p1, char p2>
double update_rules::tetra_updates<derived>::tetra(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  ffvec<p0>, ffvec<p1>, ffvec<p2>) const
{
  double p0_vec[3] = {component(p0, 0), component(p0, 1), component(p0, 2)};
  double p1_vec[3] = {component(p1, 0), component(p1, 1), component(p1, 2)};
  double p2_vec[3] = {component(p2, 0), component(p2, 1), component(p2, 2)};
  double tmp = tetra(p0_vec, p1_vec, p2_vec, u0, u1, u2, s, s0, s1, s2, h);
#if PRINT_UPDATES
  printf("tetra<%d, %d, %d>(u0 = %g, u1 = %g, u2 = %g, s = %g, "
         "s0 = %g, s1 = %g, s2 = %g, h = %g) -> %g\n",
         p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, tmp);
#endif
  return tmp;
}

template <class derived>
double update_rules::tetra_updates<derived>::tetra(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  double u[3] = {u0, u1, u2};
  double s_hat = s;
  double s_[3] = {s0, s1, s2};
  double p[3][3];
  memcpy((void *) p[0], (void *) p0, 3*sizeof(double));
  memcpy((void *) p[1], (void *) p1, 3*sizeof(double));
  memcpy((void *) p[2], (void *) p2, 3*sizeof(double));
  using cost_func_t = typename derived::cost_func;
  cost_func_t func {h, static_cast<derived const *>(this)->theta()};
  func.set_args(u, s_hat, s_, p);
  double lambda[2], F0;
  bool error;
  numopt::sqp_baryplex<cost_func_t, 3, 2> sqp;
  sqp(func, lambda, &error);
  assert(!error);
  func.set_lambda(lambda); // TODO: maybe unnecessary
  func.eval(F0);
  return F0;
}

#endif // __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__
