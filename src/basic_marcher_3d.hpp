#ifndef __BASIC_MARCHER_3D_HPP__
#define __BASIC_MARCHER_3D_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"

struct basic_marcher_3d: public marcher_3d<basic_marcher_3d, node_3d, 6>
{
  using marcher_3d::marcher_3d;
  static constexpr int num_neighbors = 6;

  double s_hat;
  node_3d * nb[num_neighbors];

EIKONAL_PRIVATE:
  virtual void update_impl(node_3d * n, node_3d ** nb, int parent, double & T);
};

#endif // __BASIC_MARCHER_3D_HPP__
