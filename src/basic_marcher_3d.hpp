#ifndef __BASIC_MARCHER_3D_HPP__
#define __BASIC_MARCHER_3D_HPP__

#include "marcher_3d.hpp"
#include "updates.line.hpp"

struct basic_marcher_3d: public marcher_3d<basic_marcher_3d, 6>,
                         public updates::line<RHR>
{
  static constexpr int num_nb = 6;
  static constexpr cost_func F_ = RHR;

  using marcher_3d::marcher_3d;

  double s_hat;
  int nb[num_nb];

OLIM_PRIVATE:
  virtual void update_impl(int lin, int * nb, int parent, double & U);
};

#endif // __BASIC_MARCHER_3D_HPP__
