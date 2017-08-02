#ifndef __OLIM8_MP0_HPP__
#define __OLIM8_MP0_HPP__

#include "moore_marcher.hpp"
#include "node.hpp"

struct olim8_mp0: public moore_marcher<node> {
  using moore_marcher::moore_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
  double get_s_est(double s, int i0, int j0, int i1, int j1);
};

#endif // __OLIM8_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
