#ifndef __MOORE_MARCHER_HPP__
#define __MOORE_MARCHER_HPP__

#include "fast_marcher.hpp"

struct moore_marcher: public fast_marcher
{
  using fast_marcher::fast_marcher;
private:
  virtual void get_valid_neighbors(size_t i, size_t j, node ** nb);
};

#endif // __MOORE_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
