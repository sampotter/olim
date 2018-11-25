#ifndef __BASIC_MARCHER_HPP__
#define __BASIC_MARCHER_HPP__

#include "marcher.hpp"
#include "node.hpp"

struct basic_marcher: marcher<basic_marcher, node, 4>
{
  static constexpr int num_neighbors = 4;

  using marcher<basic_marcher, node, num_neighbors>::marcher;

EIKONAL_PRIVATE:
  virtual void update_impl(node * n, double & T);

  double s_hat;
  node * nb[num_neighbors];
};

#endif // __BASIC_MARCHER_HPP__
