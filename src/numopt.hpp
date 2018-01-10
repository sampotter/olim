#ifndef __NUMOPT_HPP__
#define __NUMOPT_HPP__

/**
 * Numerical optimization functions
 */

#include <armadillo>
#include <functional>

#include "common.macros.hpp"

namespace numopt {
  /**
   * Returns a sorted vector of unsigned integers containing the
   * elements in [low, high) - u.
   */
  arma::uvec setdiff(unsigned low, unsigned high, arma::uvec const & u);

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
  arma::vec qpez_schur(arma::mat const & G, arma::vec const & c,
                       arma::mat const & A);
  
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
  arma::vec qpi(arma::mat const & G, arma::vec const & c, arma::mat const & A,
                arma::vec const & b, arma::vec const * x0 = nullptr,
                bool * error = nullptr, double tol = EPS(double),
                int niters = 0);

  using field_t = std::function<double(arma::vec const &)>;
  using grad_t = std::function<arma::vec(arma::vec const &)>;
  using hess_t = std::function<arma::mat(arma::vec const &)>;

  /**
   * Solve the sequential quadratic program (SQP):
   *
   *     minimize  f(x)
   *   subject to  A*x >= b
   *
   * where df is the gradient of and d2f is the Hessian of f.
   */
  arma::vec sqp(field_t const & f, grad_t const & df, hess_t const & d2f,
                arma::mat const & A, arma::vec const & b, arma::vec const & x0,
                bool * error, double tol = EPS(double), int niters = 0);
}

#endif // __NUMOPT_HPP__
