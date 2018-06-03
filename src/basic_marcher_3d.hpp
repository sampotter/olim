#ifndef __BASIC_MARCHER_3D_HPP__
#define __BASIC_MARCHER_3D_HPP__

#include "marcher_3d.hpp"
#include "node_3d.hpp"

template <class base>
struct basic_marcher_3d: public base
{
EIKONAL_PRIVATE:
  virtual void update_impl(int i, int j, int k, int parent, double & T);
};

#include "basic_marcher_3d.impl.hpp"

#endif // __BASIC_MARCHER_3D_HPP__
