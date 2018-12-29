#ifndef __NUMOPT_HPP__
#define __NUMOPT_HPP__

/**
 * Numerical optimization functions
 */

#include <src/config.hpp>

#include <functional>

#include "common.hpp"
#include "cost_funcs.hpp"

enum class line_search {BACKTRACK, HYBRID};

template <int d, int m>
void qpe_bary(double const * G, double const * c, double * x);

template <int d>
void qpi_bary(double const * G, double const * c, double const * x0,
              double * x, bool * error = nullptr,
              double tol = eps<double>, int niters = 0);

template <class cost_functor, int n, int d, line_search L = line_search::HYBRID>
struct sqp_bary {};

template <class cost_functor, line_search L>
struct sqp_bary<cost_functor, 3, 2, L> {
  void operator()(cost_functor & func, double const * xinit,
                  double * x, double * f, bool * error,
                  double tol = eps<double>, int niters = 0);
};

#include "numopt.impl.hpp"

#endif // __NUMOPT_HPP__
