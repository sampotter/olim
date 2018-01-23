#ifndef __NUMOPT_HPP__
#define __NUMOPT_HPP__

/**
 * Numerical optimization functions
 */

#include <src/config.hpp>

#if USE_ARMADILLO
#    include <armadillo>
#endif
#include <functional>

#include "common.macros.hpp"
#include "cost_funcs.hpp"

namespace numopt {
  /**
   * Returns a sorted vector of unsigned integers containing the
   * elements in [low, high) - u.
   */
#if USE_ARMADILLO
  arma::uvec setdiff(unsigned low, unsigned high, arma::uvec const & u);
#endif // USE_ARMADILLO

  /**
   * Solve the equality-constrained quadratic program:
   *
   *     minimize  x'*G*x/2 + x'*c
   *   subject to  A*x = 0
   *
   * using the Schur-complement method (see e.g. Nocedal & Wright or
   * Bertsekas. The load vector for the equality constraint is
   * zero. Requires G positive definite.
   */
#if USE_ARMADILLO
  arma::vec qpez_schur(arma::mat const & G, arma::vec const & c,
                       arma::mat const & A);
#endif // USE_ARMADILLO
  
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

  /**
   * Solve the inequality-constrained quadratic program:
   *
   *     minimize  x'*G*x/2 + x'*c
   *   subject to  A*x >= b
   *
   * using the active set method (see e.g. Nocedal & Wright or
   * Bertsekas). Equality-constrained quadratic subproblems are solved
   * using `qpe'. Requires G positive definite.
   */
#if USE_ARMADILLO
  arma::vec qpi(arma::mat const & G, arma::vec const & c, arma::mat const & A,
                arma::vec const & b, arma::vec const * x0 = nullptr,
                bool * error = nullptr, double tol = EPS(double),
                int niters = 0);
#endif // USE_ARMADILLO

  template <int d>
  void qpi_baryplex(double const * G, double const * c, double const * x0,
                    double * x, bool * error = nullptr,
                    double tol = EPS(double), int niters = 0);

#if USE_ARMADILLO
  using field_t = std::function<double(arma::vec const &)>;
  using grad_t = std::function<arma::vec(arma::vec const &)>;
  using hess_t = std::function<arma::mat(arma::vec const &)>;
#endif // USE_ARMADILLO

  /**
   * Solve the sequential quadratic program (SQP):
   *
   *     minimize  f(x)
   *   subject to  A*x >= b
   *
   * where df is the gradient of and d2f is the Hessian of f.
   */
#if USE_ARMADILLO
  arma::vec sqp(field_t const & f, grad_t const & df, hess_t const & d2f,
                arma::mat const & A, arma::vec const & b, arma::vec const & x0,
                bool * error, double tol = EPS(double), int niters = 0);
#endif // USE_ARMADILLO

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
