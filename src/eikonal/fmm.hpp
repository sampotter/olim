#pragma once

#include "line.hpp"
#include "marcher.hpp"

template <int n, ordering ord>
struct fmm:
  public eikonal::marcher<fmm<n, ord>, n, 2*n, ord>,
  public eikonal::line<RHR>
{
  static constexpr cost_func F_ = RHR;
  static constexpr int num_nb = 2*n;

  using eikonal::marcher<fmm<n, ord>, n, 2*n, ord>::marcher;

OLIM_PRIVATE:
  void update_impl(int lin, int const * nb, int parent, double & U);
};

#include "fmm.impl.hpp"
