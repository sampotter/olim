#ifndef __NUMOPT_HPP__
#define __NUMOPT_HPP__

/**
 * Numerical optimization functions
 */

#include <src/config.hpp>

#include <functional>

#include "common.macros.hpp"
#include "cost_funcs.hpp"

namespace numopt {
  /**
   * Solve the equality-constrained quadratic program:
   *
   *     minimize x'*G*x/2 + x'*c
   *   subject to a'*x = b
   *
   * where the equality constraint corresponds to the mth active
   * constraint for the barycentric simplex in R^d (of d + 1 possible
   * constraints).
   */
  template <int d, int m>
  void qpe_baryplex(double const * G, double const * c, double * x);

  template <int d>
  void qpi_baryplex(double const * G, double const * c, double const * x0,
                    double * x, bool * error = nullptr,
                    double tol = EPS(double), int niters = 0);

  template <class cost_func_t, int n, int d>
  struct sqp_baryplex {};

  template <class cost_func_t>
  struct sqp_baryplex<cost_func_t, 3, 2> {
    void operator()(cost_func_t & func, double * x, bool * error,
                    double tol = EPS(double), int niters = 0);
  };
}

#include "numopt.impl.hpp"

#endif // __NUMOPT_HPP__
