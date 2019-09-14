#pragma once

#include "line.hpp"
#include "marcher.hpp"

namespace eikonal {

template <int n, ordering ord>
struct fmm:
  public marcher<fmm<n, ord>, n, 2*n, ord>,
  public line<RHR>
{
  static constexpr cost_func F_ = RHR;
  static constexpr int num_nb = 2*n;

  using eikonal::marcher<fmm<n, ord>, n, 2*n, ord>::marcher;

  void update_impl(int lin, int const * nb, int parent, double & U);
};

template <ordering ord = ordering::COLUMN_MAJOR>
using fmm2 = fmm<2, ord>;

template <ordering ord = ordering::COLUMN_MAJOR>
using fmm3 = fmm<3, ord>;

}

#include "fmm.impl.hpp"
