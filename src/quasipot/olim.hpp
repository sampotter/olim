#pragma once

#include "../vec.hpp"
#include "marcher.hpp"

namespace quasipot {

template <ordering ord>
struct olim: public marcher<olim<ord>, 2, ord>
{
  using marcher_t = marcher<olim<ord>, 2, ord>;

  friend marcher_t;

  using marcher<olim<ord>, 2, ord>::marcher;

OLIM_PRIVATE:

  void update_impl(int lin, double & U);
};

}
