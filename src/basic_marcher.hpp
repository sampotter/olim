#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "marcher.hpp"
#include "updates.line.hpp"

struct basic_marcher: public marcher<basic_marcher, 4>,
                      public updates::line<RHR>
{
  static constexpr cost_func F_ = RHR;
  static constexpr int num_nb = 4;

  using marcher<basic_marcher, num_nb>::marcher;

  double s_hat;
  int nb[num_nb];

OLIM_PRIVATE:
  virtual void update_impl(int lin, double & U);
};

#endif // __BASIC_MARCHER_HPP__
