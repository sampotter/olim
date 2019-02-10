#pragma once

#include "marcher.hpp"
#include "updates.line.hpp"

struct basic_marcher: public marcher<basic_marcher, 2, 4>,
                      public updates::line<RHR>
{
  static constexpr cost_func F_ = RHR;
  static constexpr int num_nb = 4;

  using marcher<basic_marcher, 2, num_nb>::marcher;

OLIM_PRIVATE:
  virtual void update_impl(int lin, int const * nb, int parent, double & U);
};
