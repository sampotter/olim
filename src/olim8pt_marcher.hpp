#ifndef __OLIM8PT_MARCHER_HPP__
#define __OLIM8PT_MARCHER_HPP__

#include "moore_marcher.hpp"

struct olim8pt_marcher: public moore_marcher
{
  using moore_marcher::moore_marcher;
private:
  virtual void update_node_value_impl(size_t i, size_t j);
};

#endif // __OLIM8PT_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
