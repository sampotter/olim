#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "fast_marcher.hpp"

struct basic_marcher: public fast_marcher
{
  using fast_marcher::fast_marcher;
private:
  virtual void update_node_value_impl(size_t i, size_t j);
  virtual void update_neighbors_impl(size_t i, size_t j);
  virtual void get_valid_neighbors(size_t i, size_t j, node ** nb);
};

#endif // __BASIC_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
