#pragma once

/**
 * Numerical optimization functions
 */

#include "common.hpp"
#include "cost_funcs.hpp"
#include "vec.hpp"

enum class line_search {BACKTRACK, HYBRID};

template <int d, int m>
void qpe_bary(double const * G, vec<double, d> const & c, vec<double, d> & x);

template <int d>
void qpi_bary(double const * G, vec<double, d> const & c,
              vec<double, d> const & x0, vec<double, d> & x,
              bool * error = nullptr,
              double tol = eps<double>, int niters = 0);

template <class cost_functor, int n, int d, line_search L = line_search::HYBRID>
struct sqp_bary {};

template <class cost_functor, line_search L>
struct sqp_bary<cost_functor, 3, 2, L> {
  void operator()(cost_functor & func, vec<double, 2> const & xinit,
                  vec<double, 2> & x, double * f, bool * error,
                  double tol = eps<double>, int niters = 0);
};

#include "numopt.impl.hpp"
