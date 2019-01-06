#ifndef __UPDATES_TETRA_IMPL_HPP__
#define __UPDATES_TETRA_IMPL_HPP__

#include <src/config.hpp>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

// TODO: remove me
#include <type_traits>

#include "common.hpp"
#include "cost_funcs.hpp"
#include "numopt.hpp"

namespace updates {

template <cost_func F, int n>
inline bool should_skip(cost_functor<F, n, 2> & func, info<2> const & info) {
  if (!info.on_boundary()) {
    return false;
  }

  func.set_lambda(info.lambda);

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

// TODO: we want to cache the geometry of each workspace and pass that
// in instead of p0, p1, and p2
template <cost_func F, int n>
void
updates::tetra<F, n>::operator()(
  cost_functor<F, n, 2> & func, info<2> & info) const
{
  // TODO: turn this back on
  // if (should_skip(func, info)) {
  //   return;
  // }

  if (F == MP1) {
    bool error;
    sqp_bary<decltype(func), n, 2>()(
      func,
      nullptr, // TODO: warm start
      info.lambda,
      &info.value,
      &error);
    assert(!error);
  } else {
    auto & w = func.w;
    auto & qr = func.qr;

    // Compute A = inv(R')*du/sh here.
    double A[2] = {w.du[0]/w.sh_lam, w.du[1]/w.sh_lam};
    A[1] -= qr->r[1]*A[0]/qr->r[0];
    A[1] /= qr->r[2];
    A[0] /= qr->r[0];

    double A_dot_A = dot<2>(A, A);
    if (A_dot_A < 1) {
      double lopt = sqrt(qr->numer/(1 - A_dot_A));

      auto & lam = info.lambda;
      lam[0] = qr->Qt_p0[0] + lopt*A[0];
      lam[1] = qr->Qt_p0[1] + lopt*A[1];
      lam[1] /= -qr->r[2];
      lam[0] += qr->r[1]*lam[1];
      lam[0] /= -qr->r[0];

      if (lam[0] >= 0 && lam[1] >= 0 && lam[0] + lam[1] <= 1) {
        info.value = w.u0 + w.du[0]*lam[0] + w.du[1]*lam[1] + w.sh_lam*lopt;
      }
    }
  }
}

template <cost_func F, int n>
void
updates::tetra<F, n>::operator()(
  cost_functor_fac<F, n, 2> & func, info<2> & info) const
{
  bool error;
  sqp_bary<decltype(func), n, 2, line_search::BACKTRACK>()(
    func, nullptr, info.lambda, &info.value, &error);
  assert(!error);
}

#define __r11 bitops::R<p0, p1, p2, 0>(bitops::dim<3> {})
#define __r12 bitops::R<p0, p1, p2, 1>(bitops::dim<3> {})
#define __r22 bitops::R<p0, p1, p2, 2>(bitops::dim<3> {})
#define __numer bitops::exact_soln_numer<p0, p1, p2>(bitops::dim<3> {})
#define __Qt_p0(j) bitops::Qt_dot_p0<p0, p1, p2, j>(bitops::dim<3> {})

template <cost_func F, int n, int p0, int p1, int p2>
void
updates::tetra_bv<F, n, p0, p1, p2>::operator()(
  cost_functor_bv<F, n, p0, p1, p2> & func, info<2> & info) const
{
  if (F == MP1) {
    bool error;
    sqp_bary<decltype(func), n, 2>()(
      func,
      info.on_boundary() ? nullptr : info.lambda, // TODO: ???
      info.lambda,
      &info.value,
      &error);
    assert(!error);
  } else {
    auto & w = func.w;

    // Compute A = inv(R')*du/sh here.
    double A[2] = {w.du[0]/w.sh_lam, w.du[1]/w.sh_lam};
    A[1] -= __r12*A[0]/__r11;
    A[1] /= __r22;
    A[0] /= __r11;

    double A_dot_A = dot<2>(A, A);
    if (A_dot_A < 1) {
      double lopt = sqrt(__numer/(1 - A_dot_A));

      auto & lam = info.lambda;
      lam[0] = __Qt_p0(0) + lopt*A[0];
      lam[1] = __Qt_p0(1) + lopt*A[1];
      lam[1] /= -__r22;
      lam[0] += __r12*lam[1];
      lam[0] /= -__r11;

      if (lam[0] >= 0 && lam[1] >= 0 && lam[0] + lam[1] <= 1) {
        info.value = w.u0 + w.du[0]*lam[0] + w.du[1]*lam[1] + w.sh_lam*lopt;
      }
    }
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
