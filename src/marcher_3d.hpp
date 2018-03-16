#ifndef __MARCHER_3D_HPP__
#define __MARCHER_3D_HPP__

#include <functional>

#include "abstract_marcher.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

template <class Node>
struct marcher_3d: public abstract_marcher {
  using float_type = double;
  using node_type = Node;

  static constexpr int ndim = 3;

  marcher_3d(int height, int width, int depth, double h,
             no_speed_func_t const &);
  marcher_3d(int height, int width, int depth, double h = 1,
             std::function<double(double, double, double)> speed =
             static_cast<speed_func_3d>(default_speed_func),
             double x0 = 0.0, double y0 = 0.0, double z0 = 0.0);
  marcher_3d(int height, int width, int depth, double h,
             double const * s_cache);
  virtual ~marcher_3d();

  void add_boundary_node(int i, int j, int k, double value = 0.0);
  void add_boundary_node(double x, double y, double z, double value = 0.0);
  void add_boundary_nodes(Node const * nodes, int num_nodes);
  Node * get_node_pointer() const { return _nodes; }
  double get_speed(int i, int j, int k) const;
  double get_value(int i, int j, int k) const;
  int get_height() const { return _height; }
  int get_width() const { return _width; }
  int get_depth() const { return _depth; }
  void * get_s_cache_data() { return (void *) _s_cache; }

  Node & operator()(int i, int j, int k);
  Node const & operator()(int i, int j, int k) const;

EIKONAL_PROTECTED:
  void update(int i, int j, int k);
  void stage(int i, int j, int k);
  bool in_bounds(int i, int j, int k) const;
  bool is_valid(int i, int j, int k) const;
  double get_h() const { return _h; }

  virtual void get_valid_neighbors(int i, int j, int k, abstract_node ** nb) = 0;
  virtual void update_impl(int i, int j, int k, double & T) = 0;
  
EIKONAL_PRIVATE:
  void init();

  Node * _nodes;
  double const * _s_cache {nullptr};
  double _h {1};
  int _height, _width, _depth;
};

#include "marcher_3d.impl.hpp"

#endif // __MARCHER_3D_HPP__
