#ifndef __NUMOPT_IMPL_HPP__
#define __NUMOPT_IMPL_HPP__

#include <algorithm>
#include <cmath>

#include "hybrid.hpp"

template <int d>
inline double compute_lambda_min_symmetric(double const * A);

template <>
inline double compute_lambda_min_symmetric<2>(double const * A) {
  double const half_tr = (A[0] + A[2])/2;
  double const det = A[0]*A[2] - A[1]*A[1];
  return half_tr - sqrt(half_tr*half_tr - det);
}

template <class cost_functor>
void
sqp_bary<cost_functor, 3, 2>::operator()(
  cost_functor & func, double * x, double * f, bool * error, double tol,
  int niters)
{
  using std::max;

  if (error) *error = false;

  int k = 0;
  double x0[2], x1[2], f0, f1, lambda_min;
  double d2f[3], df[2], h[2], c[2], alpha;
  bool qpi_error;
  hybrid_status status;

  x1[0] = x1[1] = 1./3;
  func.set_lambda(x1);
  func.eval(f1);

  while (true) {
    // Compute Hessian
    func.hess(d2f);

    // ... perturb it if it isn't positive definite
    lambda_min = compute_lambda_min_symmetric<2>(d2f);
    if (lambda_min < 0) {
      d2f[0] -= 1.1*lambda_min;
      d2f[2] -= 1.1*lambda_min;
    }

    func.grad(df);

    c[0] = df[0] - d2f[0]*x1[0] - d2f[1]*x1[1];
    c[1] = df[1] - d2f[1]*x1[0] - d2f[2]*x1[1];

    x0[0] = x1[0];
    x0[1] = x1[1];

    qpi_bary<2>(d2f, c, x0, x1, &qpi_error, tol, 10);
    assert(!qpi_error);

    // Compute descent step: h = x1 - x0 = x_{n+1} - x_{n}
    sub<2>(x1, x0, h);

    // Apply Wilkinson hybrid method to minimize f(x0 + alpha*h) for
    // alpha st 0 <= alpha <= 1. We do this by finding the zero of
    // df(x0 + alpha*h)'*alpha over the same interval.
    auto const step = [&] (double alpha) {
      axpy<2>(alpha, h, x0, x1);
      func.set_lambda(x1);
      func.grad(df);
      return dot<2>(df, h);
    };
    std::tie(alpha, status) = hybrid(step, 0., 1., tol);
    if (status == hybrid_status::DEGENERATE) {
      alpha = 1;
    }

    axpy<2>(alpha, h, x0, x1);
    f0 = f1;
    func.set_lambda(x1);
    func.eval(f1);

    if (fabs(f1 - f0) <= tol*fmax(f0, f1) + tol) {
      break;
    }

    if (fmax(fabs(x1[0] - x0[0]), fabs(x1[1] - x0[1])) <=
        tol*fmax(fmax(x0[0], x0[1]), fmax(x1[0], x1[1])) + tol) {
      break;
    }

    // Check if we've reached our max number of iterations
    if (++k == niters) {
      if (error) *error = true;
      break;
    }
  }

  if (x != nullptr) {
    x[0] = x1[0];
    x[1] = x1[1];
  }

  if (f != nullptr) {
    *f = f1;
  }
}

#undef __compute_lambda_min

#endif // __NUMOPT_IMPL_HPP__
