#ifndef __BASIC_MARCHER_3D_HPP__
#define __BASIC_MARCHER_3D_HPP__

#include "neumann_marcher_3d.hpp"
#include "node_3d.hpp"

struct basic_marcher_3d: public neumann_marcher_3d<node_3d>
{
  using neumann_marcher_3d::neumann_marcher_3d;
EIKONAL_PRIVATE:
  virtual void update_impl(int i, int j, int k, double & T);
};

#endif // __BASIC_MARCHER_3D_HPP__
