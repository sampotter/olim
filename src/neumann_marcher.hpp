#ifndef __NEUMANN_MARCHER_HPP__
#define __NEUMANN_MARCHER_HPP__

#include "fast_marcher.hpp"

struct neumann_marcher: public fast_marcher {
  using fast_marcher::fast_marcher;
protected:
  virtual void get_valid_neighbors(size_t i, size_t j, node** nb);
private:
  virtual void stage_neighbors_impl(size_t i, size_t j);
};

#endif // __NEUMANN_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
