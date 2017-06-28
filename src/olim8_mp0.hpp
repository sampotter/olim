#ifndef __OLIM8_MP0_HPP__
#define __OLIM8_MP0_HPP__

#include "moore_marcher.hpp"

struct olim8_mp0: public moore_marcher
{
  using moore_marcher::moore_marcher;
private:
  virtual void update_node_value_impl(size_t i, size_t j, double & T);
  double get_s_est(double s, size_t i0, size_t j0, size_t i1, size_t j1);
};

#endif // __OLIM8_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
