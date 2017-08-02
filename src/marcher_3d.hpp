#ifndef __MARCHER_3D_HPP__
#define __MARCHER_3D_HPP__

#include "abstract_marcher.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

template <class Node>
struct marcher_3d: public abstract_marcher
{
  void add_boundary_node(int i, int j, int k, double value = 0.0);
  double get_value(int i, int j, int k) const;

protected:
  Node & operator()(int i, int j, int k);
  Node const & operator()(int i, int j, int k) const;
  void update(int i, int j, int k);
  void stage(int i, int j, int k);
  bool in_bounds(int i, int j, int k) const;
  bool is_valid(int i, int j, int k) const;
  double S(int i, int j, int k);
  double get_h() const { return _h; }

  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb) = 0;
  
private:
  virtual void update_node_value_impl(int i, int j, int k, double & T) = 0;
  virtual void stage_neighbors_impl(abstract_node * n) = 0;

  Node * _nodes;
  speed_func_3d _S {default_speed_func_3d};
  double _h {1};
  double _x0 {0}, _y0 {0}, _z0 {0};
  int _height, _width, _depth;
};

#include "marcher_3d.impl.hpp"

#endif // __MARCHER_3D_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
