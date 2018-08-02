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
#include <type_traits>

#include "common.defs.hpp"
#include "common.macros.hpp"
#include "numopt.hpp"
#include "update_rules.tetra_updates.util.hpp"

#define __max3(a, b, c) std::max(a, std::max(b, c))

#define __theta() static_cast<derived const *>(this)->theta()

template <class derived>
template <char p0, char p1, char p2>
update_info<2> update_rules::tetra_updates<derived>::tetra(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  ffvec<p0>, ffvec<p1>, ffvec<p2>) const
{
  double u[3] = {u0, u1, u2};
  double s_hat = s;
  double s_[3] = {s0, s1, s2};

  using cost_func_t = typename derived::template cost_func_bv<p0, p1, p2>;

  cost_func_t func {h, __theta()};
  func.set_args(u, s_hat, s_);

  update_info<2> update;
  bool error;
  numopt::sqp_baryplex<cost_func_t, 3, 2> sqp;
  sqp(func, update.lambda, &error);
  assert(!error);

  // TODO: awful hack for now---need to fix the way we've organized
  // the cost functions
  if (std::is_same<derived, mp0_tetra_updates>::value) {
    F1_bv<p0, p1, p2, 2> eval_func {h, __theta()};
    eval_func.set_args(u, s_hat, s_);
    eval_func.set_lambda(update.lambda);
    eval_func.eval(update.value);
  } else {
    func.set_lambda(update.lambda); // TODO: maybe unnecessary
    func.eval(update.value);
  }

#if PRINT_UPDATES
  printf("tetra<%d, %d, %d>(u0 = %g, u1 = %g, u2 = %g, s = %g, "
         "s0 = %g, s1 = %g, s2 = %g, h = %g) -> %g\n",
         p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, value);
#endif

  assert(update.value + EPS(double) >= __max3(u0, u1, u2));

  return update;
}

template <class derived>
update_info<2> update_rules::tetra_updates<derived>::tetra(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  double const * p0, double const * p1, double const * p2,
  double const * p_fac, double s_fac) const  
{
  using cost_func_t = typename derived::template cost_func<3, 2>;

  using eval_func_t = typename std::conditional<
    std::is_same<derived, mp0_tetra_updates>::value,
    F1<3, 2>,
    cost_func_t
  >::type;

  using factored_cost_func_t =
    typename derived::template factored_cost_func<3, 2>;

  // TODO: we want to avoid having to pack everything into these
  // arrays, this is a waste
  double u[3] = {u0, u1, u2};
  double s_hat = s;
  double s_[3] = {s0, s1, s2};
  double p[3][3];

  // TODO: we want to avoid doing this, this is totally unnecessary
  memcpy((void *) p[0], (void *) p0, 3*sizeof(double));
  memcpy((void *) p[1], (void *) p1, 3*sizeof(double));
  memcpy((void *) p[2], (void *) p2, 3*sizeof(double));
  
  factored_cost_func_t func {h, __theta()};
  func.set_args(u, s_hat, s_, p, p_fac, s_fac);

  update_info<2> update;
  bool error;
  numopt::sqp_baryplex<factored_cost_func_t, 3, 2> sqp;
  sqp(func, update.lambda, &error);
  assert(!error);

  eval_func_t eval_func {h, __theta()};
  eval_func.set_args(u, s_hat, s_, p);
  eval_func.set_lambda(update.lambda);
  eval_func.eval(update.value);

#if PRINT_UPDATES
#  error Not implemented yet
#endif

  assert(update.value >= __max3(u0, u1, u2));

  return update;
}

template <class derived>
update_info<2> update_rules::tetra_updates<derived>::tetra(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  // TODO: we want to avoid having to pack everything into these
  // arrays, this is a waste
  double u[3] = {u0, u1, u2};
  double s_hat = s;
  double s_[3] = {s0, s1, s2};
  double p[3][3];

  // TODO: we want to avoid doing this, this is totally unnecessary
  memcpy((void *) p[0], (void *) p0, 3*sizeof(double));
  memcpy((void *) p[1], (void *) p1, 3*sizeof(double));
  memcpy((void *) p[2], (void *) p2, 3*sizeof(double));

  using cost_func_t = typename derived::template cost_func<3, 2>;

  cost_func_t func {h, __theta()};
  func.set_args(u, s_hat, s_, p);

  update_info<2> update;
  bool error;
  numopt::sqp_baryplex<cost_func_t, 3, 2> sqp;
  sqp(func, update.lambda, &error);
  assert(!error);

  // TODO: awful hack for now---need to fix the way we've organized
  // the cost functions
  if (std::is_same<derived, mp0_tetra_updates>::value) {
    F1<3, 2> eval_func {h, __theta()};
    eval_func.set_args(u, s_hat, s_, p);
    eval_func.set_lambda(update.lambda);
    eval_func.eval(update.value);
  } else {
    func.set_lambda(update.lambda); // TODO: maybe unnecessary
    func.eval(update.value);
  }

#if PRINT_UPDATES
  printf("tetra(u0 = %g, u1 = %g, u2 = %g, s = %g, "
         "s0 = %g, s1 = %g, s2 = %g, h = %g) -> %g\n",
         u0, u1, u2, s, s0, s1, s2, h, value);
#endif

  assert(update.value >= __max3(u0, u1, u2));

  return update;
}

#undef __theta

#endif // __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__
