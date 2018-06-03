#ifndef __BASIC_MARCHER_3D_HPP__
#define __BASIC_MARCHER_3D_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"

struct basic_marcher_3d: public marcher_3d<basic_marcher_3d, node_3d> {
  using marcher_3d::marcher_3d;
  static constexpr int nneib = 6;
EIKONAL_PRIVATE:
  virtual void update_impl(int i, int j, int k, int parent,
                           abstract_node ** nb, double & T);
};

#endif // __BASIC_MARCHER_3D_HPP__
