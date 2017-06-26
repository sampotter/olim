#ifndef __OLIM_8PT_MP0_HPP__
#define __OLIM_8PT_MP0_HPP__

#include "moore_marcher.hpp"

struct olim_8pt_mp0: public moore_marcher
{
  using moore_marcher::moore_marcher;
private:
  virtual void update_node_value_impl(size_t i, size_t j, double & T);
  double get_s_est(double s, double x0, double y0, double x1, double y1);
};

#endif // __OLIM_8PT_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
