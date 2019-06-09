#pragma once

#include "line.hpp"
#include "marcher.hpp"

template <int N>
struct fmm: public marcher<fmm<N>, N, 2*N>, public eikonal::line<RHR>
{
  static constexpr cost_func F_ = RHR;
  static constexpr int num_nb = 2*N;

  using marcher<fmm<N>, N, 2*N>::marcher;

OLIM_PRIVATE:
  virtual void update_impl(int lin, int const * nb, int parent, double & U);
};

#include "fmm.impl.hpp"