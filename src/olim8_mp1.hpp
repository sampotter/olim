#ifndef __OLIM8_MP1_HPP__
#define __OLIM8_MP1_HPP__

#include "moore_marcher.hpp"

struct olim8_mp1: public moore_marcher
{
  using moore_marcher::moore_marcher;
private:
  virtual void update_node_value_impl(size_t i, size_t j, double & T);
};

#endif // __OLIM8_RHR_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
