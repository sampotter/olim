#ifndef __NUMOPT_HPP__
#define __NUMOPT_HPP__

/**
 * Numerical optimization functions
 */

#include <armadillo>

#include "common.macros.hpp"

namespace numopt {
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
  arma::vec qpez(arma::mat const & G, arma::vec const & c, arma::mat const & A);
  
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
                arma::vec const & b, double tol = EPS(double), int niters = 0);
}

#include "numopt.impl.hpp"

#endif // __NUMOPT_HPP__
