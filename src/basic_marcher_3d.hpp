#pragma once

#include "marcher.hpp"
#include "updates.line.hpp"

struct basic_marcher_3d: public marcher<basic_marcher_3d, 3, 6>,
                         public updates::line<RHR>
{
  static constexpr int num_nb = 6;
  static constexpr cost_func F_ = RHR;

  using marcher::marcher;

OLIM_PRIVATE:
  virtual void update_impl(int lin, int const * nb, int parent, double & U);
};
