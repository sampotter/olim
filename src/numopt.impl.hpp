#ifndef __NUMOPT_IMPL_HPP__
#define __NUMOPT_IMPL_HPP__

#include <math.h>

#include "hybrid.hpp"

template <int d>
inline double compute_lambda_min_symmetric(double const * A);

template <>
inline double compute_lambda_min_symmetric<2>(double const * A) {
  double const half_tr = (A[0] + A[2])/2;
  double const det = A[0]*A[2] - A[1]*A[1];
  return half_tr - sqrt(half_tr*half_tr - det);
}

template <class cost_functor, line_search L>
void
sqp_bary<cost_functor, 3, 2, L>::operator()(
  cost_functor & func, vec<double, 2> const & xinit,
  vec<double, 2> & x, double * f, bool * error, double tol, int niters)
{
  if (error) *error = false;

  int k = 0;
  double f0, f1, lambda_min;
  double d2f[3], alpha;
  vec<double, 2> x0, x1, df, h, c;
  bool qpi_error;
  hybrid_status status;

  x1 = xinit;

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
    h = x1 - x0;

    if (L == line_search::BACKTRACK) {
      alpha = 1;
      double lhs, rhs = f1 + 1e-4*alpha*(df*h);
    repeat:
      x1 = x0 + alpha*h;
      func.set_lambda(x1);
      func.eval(lhs);
      if (alpha > tol && lhs > rhs) {
        alpha /= 2;
        goto repeat;
      }
    } else if (L == line_search::HYBRID) {
      // Apply Wilkinson's hybrid method to minimize f(x0 + alpha*h) for
      // alpha st 0 <= alpha <= 1. We do this by finding the zero of
      // df(x0 + alpha*h)'*alpha over the same interval.
      auto const step = [&] (double alpha) {
        x1 = x0 + alpha*h;
        func.set_lambda(x1);
        func.grad(df);
        return df*h;
      };
      std::tie(alpha, status) = hybrid(step, 0., 1., tol);
      if (status == hybrid_status::DEGENERATE) {
        alpha = 1;
      }
    }

    x1 = x0 + alpha*h;
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

  x = x1;

  if (f != nullptr) {
    *f = f1;
  }
}

#undef __compute_lambda_min

#endif // __NUMOPT_IMPL_HPP__
