#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "marcher.hpp"
#include "node.hpp"
#include "updates.line.hpp"

struct basic_marcher: public marcher<basic_marcher, node, 4>,
                      public updates::line<RHR>
{
  static constexpr cost_func F_ = RHR;
  static constexpr int num_neighbors = 4;

  using marcher<basic_marcher, node, num_neighbors>::marcher;

  double s_hat;
  node * nb[num_neighbors];

EIKONAL_PRIVATE:
  virtual void update_impl(node * n, double & T);
};

#endif // __BASIC_MARCHER_HPP__
