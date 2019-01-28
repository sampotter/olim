#pragma once

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

// TODO: this should (probably) eventually become an olim3d_hu member
// function, since that's the only class that uses it...
// TODO: simplify interface using template template parameter...?
template <cost_func F, int n>
inline bool should_skip(cost_functor<F, n, 2> & func, info<2> const & info) {
  // TODO: instead, `assert(info.on_boundary())' here, since if
  // in_interior() is true, then this is false... so a bit of a waste!
  if (info.in_interior()) {
    return false;
  }

  func.set_lambda(info.lambda);

  vec<double, 2> df;
  double d2f[3];
  func.grad(df);
  func.hess(d2f);

  vec<double, 2> mu;
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
  if (F == MP1) {
    bool error;
    sqp_bary<decltype(func), n, 2>()(
      func,
      info.lambda,
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
    func, info.lambda, info.lambda, &info.value, &error);
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
      info.lambda,
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
