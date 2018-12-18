#ifndef __NUMOPT_HPP__
#define __NUMOPT_HPP__

/**
 * Numerical optimization functions
 */

#include <src/config.hpp>

#include <functional>

#include "common.macros.hpp"
#include "cost_funcs.hpp"

template <int d, int m>
void qpe_bary(double const * G, double const * c, double * x);

template <int d>
void qpi_bary(double const * G, double const * c, double const * x0,
              double * x, bool * error = nullptr,
              double tol = EPS(double), int niters = 0);

template <class cost_functor, int n, int d>
struct sqp_bary {};

template <class cost_functor>
struct sqp_bary<cost_functor, 3, 2> {
  void operator()(cost_functor & func, double const * xinit,
                  double * x, double * f, bool * error,
                  double tol = EPS(double), int niters = 0);
};

#include "numopt.impl.hpp"

#endif // __NUMOPT_HPP__
