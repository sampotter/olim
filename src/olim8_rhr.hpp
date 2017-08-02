#ifndef __OLIM8_RHR_HPP__
#define __OLIM8_RHR_HPP__

#include "moore_marcher.hpp"
#include "node.hpp"

struct olim8_rhr: public moore_marcher<node>
{
  using moore_marcher::moore_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
  double get_alpha(double u0, double u1, double sh) const;
};

#endif // __OLIM8_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
