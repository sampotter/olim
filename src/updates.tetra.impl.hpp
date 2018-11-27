#ifndef __UPDATES_TETRA_IMPL_HPP__
#define __UPDATES_TETRA_IMPL_HPP__

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
#include "cost_funcs.hpp"
#include "numopt.hpp"
#include "updates.tetra.utils.hpp"

template <cost_func F, int n>
updates::info<2>
updates::tetra<F, n>::operator()(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  F_wkspc<F, 2> w;
  set_args<F, n>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor<F, 2> func {w, p0, p1, p2};

  info<2> info;
  bool error;
  sqp_bary<decltype(func), n, 2>()(func, info.lambda, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  } else {
    // TODO: we're doing an unnecessary eval here: we could reorganize
    // things so that we're using the most recent eval done by
    // sqp... i.e., have sqp write over info.value internally
    // eval(func.w, info.value);
    func.eval(info.value);
  }

#if PRINT_UPDATES
  printf("tetra(u0 = %g, u1 = %g, u2 = %g, s = %g, "
         "s0 = %g, s1 = %g, s2 = %g, h = %g) -> %g\n",
         u0, u1, u2, s, s0, s1, s2, h, info.value);
#endif

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
  cost_functor_fac<F, 2> func {w, p0, p1, p2, p_fac};
  
  info<2> info;
  bool error;
  sqp_bary<decltype(func), n, 2>()(func, info.lambda, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  } else {
    func.eval(info.value);
  }

#if PRINT_UPDATES
#  error Not implemented yet
#endif

  return info;
}

template <cost_func F, int n, char p0, char p1, char p2>
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
  sqp_bary<decltype(func), n, 2>()(func, info.lambda, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  } else {
    // TODO: I think we should actually be okay to remove this---try
    // doing so once tests are stabilized after this big change
    // eval(w, info.value);
    func.eval(info.value);
  }

#if PRINT_UPDATES
  printf("tetra<%d, %d, %d>(u0 = %g, u1 = %g, u2 = %g, s = %g, "
         "s0 = %g, s1 = %g, s2 = %g, h = %g) -> %g\n",
         p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, info.value);
#endif

  return info;
}

// TODO: we aren't actually using this overload yet...
// template <cost_func F, char p0, char p1, char p2>
// updates::info<2>
// updates::tetra_bv<F, p0, p1, p2>::operator()(F_fac_wkspc<F, 2> & fw) const
// { ... }

#endif // __UPDATES_TETRA_IMPL_HPP__
