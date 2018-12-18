#ifndef __UPDATES_TETRA_IMPL_HPP__
#define __UPDATES_TETRA_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <type_traits>

#include "common.defs.hpp"
#include "common.macros.hpp"
#include "cost_funcs.hpp"
#include "numopt.hpp"

template <cost_func F, int n>
updates::info<2>
updates::tetra<F, n>::operator()(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  F_wkspc<F, 2> w;
  set_args<F, n>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor<F, 3> func {w, p0, p1, p2};

  info<2> info;
  bool error;
  sqp_bary<decltype(func), n, 2>()(func, nullptr, info.lambda, &info.value, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  }

  return info;
}

template <cost_func F, int n>
updates::info<2>
updates::tetra<F, n>::operator()(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  double const * p_fac, double s_fac) const  
{
  F_fac_wkspc<F, 2> w;
  set_args<F, n>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, p_fac, s_fac);
  cost_functor_fac<F, 3> func {w, p0, p1, p2, p_fac};
  
  info<2> info;
  bool error;
  sqp_bary<decltype(func), n, 2, line_search::BACKTRACK>()(
    func, nullptr, info.lambda, &info.value, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  }

  return info;
}

template <cost_func F, int n, int p0, int p1, int p2>
updates::info<2>
updates::tetra_bv<F, n, p0, p1, p2>::operator()(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  F_wkspc<F, 2> w;
  set_args<F, n, p0, p1, p2>(w, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor_bv<F, n, p0, p1, p2> func {w};

  info<2> info;
  bool error;
  sqp_bary<decltype(func), n, 2>()(
    func, info.is_degenerate() ? nullptr : info.lambda,
    info.lambda, &info.value, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  }

  return info;
}

// TODO: we aren't actually using this overload yet...
// template <cost_func F, int p0, int p1, int p2>
// updates::info<2>
// updates::tetra_bv<F, p0, p1, p2>::operator()(F_fac_wkspc<F, 2> & fw) const
// { ... }

#endif // __UPDATES_TETRA_IMPL_HPP__
