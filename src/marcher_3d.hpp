#ifndef __MARCHER_3D_HPP__
#define __MARCHER_3D_HPP__

#include "abstract_marcher.hpp"
#include "heap.hpp"
#include "node_3d.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

struct marcher_3d: public abstract_marcher
{
  void add_boundary_node(int i, int j, int k, double value = 0.0);
  double get_value(int i, int j, int k) const;

protected:
  node_3d & operator()(int i, int j, int k);
  node_3d const & operator()(int i, int j, int k) const;
  void update_node_value(int i, int j, int k);
  void stage_neighbor(int i, int j, int k);
  bool in_bounds(int i, int j, int k) const;
  bool is_valid(int i, int j, int k) const;
  double S(int i, int j, int k);

  virtual void get_valid_neighbors(int i, int j, int k, node_3d ** nb) = 0;
  
private:
  virtual void update_node_value_impl(int i, int j, int k, double & T) = 0;
  virtual void stage_neighbors_impl(abstract_node * n) = 0;

  node_3d * _nodes;
  heap<node_3d> _heap;
  speed_func_3d _S {default_speed_func_3d};
  double _x0 {0}, _y0 {0}, _z0 {0};
  int _height, _width, _depth;
};

#endif // __MARCHER_3D_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
