#ifndef __MARCHER_3D_HPP__
#define __MARCHER_3D_HPP__

#include "abstract_marcher.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

template <class Node>
struct marcher_3d: public abstract_marcher
{
  using node_type = Node;

  marcher_3d(
    int height, int width, int depth, double h = 1,
    std::function<double(double, double, double)> speed = default_speed_func_3d,
    double x0 = 0.0, double y0 = 0.0, double z0 = 0.0);
  marcher_3d(int height, int width, int depth, double h,
             std::unique_ptr<double[]> s_cache);

  void add_boundary_node(int i, int j, int k, double value = 0.0);
  void add_boundary_node(double x, double y, double z, double value = 0.0);
  void add_boundary_nodes(Node const * nodes, int num_nodes);
  double get_value(int i, int j, int k) const;

protected:
  Node & operator()(int i, int j, int k);
  Node const & operator()(int i, int j, int k) const;
  void update(int i, int j, int k);
  void stage(int i, int j, int k);
  bool in_bounds(int i, int j, int k) const;
  double speed(int i, int j, int k);
  bool is_valid(int i, int j, int k) const;
  double get_h() const { return _h; }

  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb) = 0;
  virtual void update_impl(int i, int j, int k, double & T) = 0;
  
private:
  void init();

  Node * _nodes;
  std::unique_ptr<double[]> _s_cache {nullptr};
  double _h {1};
  double _x0 {0}, _y0 {0}, _z0 {0};
  int _height, _width, _depth;
};

#include "marcher_3d.impl.hpp"

#endif // __MARCHER_3D_HPP__
