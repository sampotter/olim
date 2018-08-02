#ifndef __NUMOPT_IMPL_HPP__
#define __NUMOPT_IMPL_HPP__

#include <algorithm>
#include <cmath>

#include "hybrid.hpp"

#define __compute_lambda_min() do {                                 \
    double half_tr = (G[0] + G[2])/2, det = G[0]*G[2] - G[1]*G[1];  \
    lambda_min = half_tr - sqrt(half_tr*half_tr - det);             \
  } while (0)                                                       \

template <class cost_func_t>
void
numopt::sqp_baryplex<cost_func_t, 3, 2>::operator()(
  cost_func_t & func, double * x, bool * error, double tol, int niters)
{
  using std::max;

  if (error) *error = false;

  double G[3], x0[2], x1[2] = {1./3, 1./3}, c[2], g[2], f0, f1,
    lambda_min, qpi_tol, c1 = 1e-4, alpha;
  bool qpi_error, found_opt;
  int k = 0, qpi_niters = 10;
  hybrid_status status;

  (void) c1;

  func.set_lambda(x1);
  func.eval(f1);

  while (true) {
    // Compute Hessian and perturb it if it isn't positive definite
    func.hess(G);
    __compute_lambda_min();
    if (lambda_min < 0) {
      G[0] -= 1.1*lambda_min;
      G[2] -= 1.1*lambda_min;
    }

    // Compute load vector for quadratic program
    func.grad(c);
    c[0] -= G[0]*x1[0] + G[1]*x1[1];
    c[1] -= G[1]*x1[0] + G[2]*x1[1];

    // Compute descent direction by solving inequality-constrained
    // quadratic program
    found_opt = false;
    qpi_tol = tol;
    while (!found_opt) {
      qpi_baryplex<2>(G, c, x1, g, &qpi_error, qpi_tol, qpi_niters);
      if (qpi_error) qpi_tol *= 10;
      else found_opt = true;
    }
    g[0] -= x1[0];
    g[1] -= x1[1];

    auto const step_sel = [&] (double alpha) {
      double tmp[2] = {x1[0] + alpha*g[0], x1[1] + alpha*g[1]};
      func.set_lambda(tmp);
      func.grad(tmp);
      return g[0]*tmp[0] + g[1]*tmp[1];
    };
    std::tie(alpha, status) = hybrid(step_sel, 0., 1., tol);
    if (status == hybrid_status::DEGENERATE && fabs(alpha) < tol) {
      alpha = 1;
    }
    // if (alpha < 1) { printf("alpha = %g\n", alpha); }

    // Save current values for next iteration
    x0[0] = x1[0];
    x0[1] = x1[1];
    x1[0] += alpha*g[0];
    x1[1] += alpha*g[1];
    f0 = f1;
    func.set_lambda(x1);
    func.eval(f1);

    // Check for convergence
    if (max(fabs(x1[0] - x0[0]), fabs(x1[1] - x0[1]))/fmax(
          fmax(x0[0], x0[1]),
          fmax(x1[0], x1[1])) < tol ||
        fabs(f1 - f0)/fmax(f0, f1) < tol) {
      break;
    }

    // Check if we've reached our max number of iterations
    if (++k == niters) {
      if (error) *error = true;
      break;
    }
  }

  x[0] = x1[0];
  x[1] = x1[1];
}

#undef __compute_lambda_min

#endif // __NUMOPT_IMPL_HPP__
