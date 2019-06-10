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

OLIM_PRIVATE:
  void update_impl(int lin, int const * nb, int parent, double & U);
};

}

#include "fmm.impl.hpp"
