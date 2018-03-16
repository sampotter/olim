#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "neumann_marcher.hpp"
#include "node.hpp"

struct basic_marcher: public neumann_marcher<node> {
  using neumann_marcher::neumann_marcher;
EIKONAL_PRIVATE:
  virtual void update_impl(int i, int j, double & T);
};

#endif // __BASIC_MARCHER_HPP__
