#ifndef __UPDATES_TETRA_IMPL_HPP__
#define __UPDATES_TETRA_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <type_traits>

#include "common.hpp"
#include "cost_funcs.hpp"
#include "numopt.hpp"

namespace updates {

template <cost_func F>
inline bool should_skip(cost_functor<F, 3> const & func, info<2> const & info) {
  if (!info.on_boundary()) {
    return false;
  }

  double df[2], d2f[3];
  func.grad(df);
  func.hess(d2f);

  double mu[2];
  int k;
  lagmults<2>(info.lambda, df, d2f, mu, &k);

  if (k > 0 && mu[0] < 0) return false;
  if (k > 1 && mu[1] < 0) return false;
  return true;
}

}

template <cost_func F, int n>
void
updates::tetra<F, n>::operator()(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  info<2> & info) const
{
  F_wkspc<F, 2> w;
  set_args<F, n>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);

  cost_functor<F, 3> func {w};
  func.set_lambda(info.lambda);

  if (should_skip(func, info)) {
    return;
  }

  bool error;
  sqp_bary<decltype(func), n, 2>()(func, nullptr, info.lambda, &info.value, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  }
}

template <cost_func F, int n>
void
updates::tetra<F, n>::operator()(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  double const * p_fac, double s_fac,
  info<2> & info) const
{
  F_fac_wkspc<F, 2> w;
  set_args<F, n>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, p_fac, s_fac);
  cost_functor_fac<F, 3> func {w, p0, p1, p2, p_fac};
  
  bool error;
  sqp_bary<decltype(func), n, 2, line_search::BACKTRACK>()(
    func, nullptr, info.lambda, &info.value, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  }
}

#define __r11 bitops::R<p0, p1, p2, 0>(bitops::dim<3> {})
#define __r12 bitops::R<p0, p1, p2, 1>(bitops::dim<3> {})
#define __r22 bitops::R<p0, p1, p2, 2>(bitops::dim<3> {})
#define __numer bitops::exact_soln_numer<p0, p1, p2>(bitops::dim<3> {})
#define __Qt_p0(j) bitops::Qt_dot_p0<p0, p1, p2, j>(bitops::dim<3> {})

template <cost_func F, int n, int p0, int p1, int p2>
void
updates::tetra_bv<F, n, p0, p1, p2>::operator()(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  info<2> & info) const
{
  if (F == MP0 || F == RHR) {
    double sh = (F == RHR ? s : (s + (s0 + s1 + s2)/3)/2)*h;
    double du[2] = {u1 - u0, u2 - u0};

    // Compute q = inv(R')*du/sh here.
    double q[2] = {du[0]/sh, du[1]/sh};
    q[1] -= __r12*q[0]/__r11;
    q[1] /= __r22;
    q[0] /= __r11;

    double q_dot_q = dot<2>(q, q);
    if (q_dot_q < 1) {
      double lopt = sqrt(__numer/(1 - q_dot_q));

      auto & lam = info.lambda;
      lam[0] = __Qt_p0(0) + lopt*q[0];
      lam[1] = __Qt_p0(1) + lopt*q[1];
      lam[1] /= -__r22;
      lam[0] += __r12*lam[1];
      lam[0] /= -__r11;

      if (lam[0] >= 0 && lam[1] >= 0 && lam[0] + lam[1] <= 1) {
        info.value = u0 + du[0]*lam[0] + du[1]*lam[1];
        if (F == RHR) {
          info.value += lopt*sh;
        } else {
          info.value += lopt*(s + s0 + (s1 - s0)*lam[0] + (s2 - s0)*lam[1])*h/2;
        }
        return;
      }
    }
  }

  F_wkspc<F, 2> w;
  set_args<F, n, p0, p1, p2>(w, u0, u1, u2, s, s0, s1, s2, h);
  cost_functor_bv<F, n, p0, p1, p2> func {w};

  bool error;
  sqp_bary<decltype(func), n, 2>()(
    func, info.on_boundary() ? nullptr : info.lambda,
    info.lambda, &info.value, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    eval_mp1_fix(w, s, s0, s1, s2, h, info.lambda, info.value);
  }
}

#undef __r11
#undef __r12
#undef __r22
#undef __numer
#undef __Qt_p0

// TODO: we aren't actually using this overload yet...
// template <cost_func F, int p0, int p1, int p2>
// updates::info<2>
// updates::tetra_bv<F, p0, p1, p2>::operator()(F_fac_wkspc<F, 2> & fw) const
// { ... }

#endif // __UPDATES_TETRA_IMPL_HPP__
