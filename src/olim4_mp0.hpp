#ifndef __OLIM4_MP0_HPP__
#define __OLIM4_MP0_HPP__

#include "neumann_marcher.hpp"

struct olim4_mp0: public neumann_marcher
{
  using neumann_marcher::neumann_marcher;
private:
  virtual void update_node_value_impl(int i, int j, double & T);
};

#endif // __OLIM4_MP0_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End: