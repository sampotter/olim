#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "marcher.hpp"
#include "node.hpp"

struct basic_marcher: marcher<basic_marcher, node> {
  static constexpr int nneib = 4;
  using marcher::marcher;
EIKONAL_PRIVATE:
  virtual void update_impl(int i, int j, node ** nb, double & T);
};

#endif // __BASIC_MARCHER_HPP__
