#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "neumann_marcher.hpp"

struct basic_marcher: public neumann_marcher
{
  using neumann_marcher::neumann_marcher;
private:
  virtual void update_node_value_impl(size_t i, size_t j, double & T);
};

#endif // __BASIC_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
