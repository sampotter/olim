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
  if (info.in_interior()) {
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
    direct_solve<F, n>(func.w, func.qr, info.lambda, info.value);
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

template <cost_func F, int n, int p0, int p1, int p2>
void
updates::tetra_bv<F, n, p0, p1, p2>::operator()(
  cost_functor_bv<F, n, p0, p1, p2> & func, info<2> & info) const
{
  if (F == MP1) {
    bool error;
    sqp_bary<decltype(func), n, 2>()(
      func,
      info.inbounds() ? info.lambda : nullptr,
      info.lambda,
      &info.value,
      &error);
    assert(!error);
  } else {
    direct_solve<F, n, p0, p1, p2>(func.w, info.lambda, info.value);
  }
}

// TODO: we aren't actually using this overload yet...
// template <cost_func F, int p0, int p1, int p2>
// updates::info<2>
// updates::tetra_bv<F, p0, p1, p2>::operator()(F_fac_wkspc<F, 2> & fw) const
// { ... }

#endif // __UPDATES_TETRA_IMPL_HPP__
