#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "neumann_marcher.hpp"
#include "node.hpp"

struct basic_marcher: public neumann_marcher<node> {
  using neumann_marcher::neumann_marcher;
private:
  virtual void update_impl(int i, int j, double & T);
#if TRIAL_NODE_OPTIMIZATION
  virtual void update_impl(int i, int j, int src, double & T);
#endif // TRIAL_NODE_OPTIMIZATION
};

#endif // __BASIC_MARCHER_HPP__
