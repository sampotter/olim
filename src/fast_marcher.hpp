#ifndef __FAST_MARCHER_HPP__
#define __FAST_MARCHER_HPP__

#include "marcher.hpp"
#include "node.hpp"

struct fast_marcher: public marcher<node> {
  using marcher<node>::marcher;
};

#endif // __FAST_MARCHER_HPP__
