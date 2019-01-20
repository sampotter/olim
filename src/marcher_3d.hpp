#ifndef __MARCHER_3D_HPP__
#define __MARCHER_3D_HPP__

// TODO: try to remove this
#include <functional>

#include "abstract_marcher.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

template <class base, class node, int num_neighbors>
struct marcher_3d: public abstract_marcher {
  // These are for use with our pybind11 bindings. They aren't used
  // internally.
  using float_type = double;
  using node_type = node;

  static constexpr int ndim = 3;

  marcher_3d();
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
  void add_boundary_node(
    double x, double y, double z, double s, double value = 0.0);
  void add_boundary_nodes(node const * nodes, int num_nodes);
  void add_boundary_nodes(node const * const * nodes, int num_nodes);
  void set_node_fac_center(
    int i, int j, int k, typename node::fac_center const * fac);

  node * get_node_pointer() const { return _nodes; }
  double get_speed(int i, int j, int k) const;
  double get_value(int i, int j, int k) const;
  int get_height() const { return _height; }
  int get_width() const { return _width; }
  int get_depth() const { return _depth; }
  void * get_s_cache_data() { return (void *) _s_cache; }

  node & operator()(int i, int j, int k);
  node const & operator()(int i, int j, int k) const;

EIKONAL_PROTECTED:
  inline int linear_index(int i, int j, int k) const {
    return _height*(_width*k + j) + i; // column-major
  }

  bool in_bounds(int i, int j, int k) const;
  bool is_valid(int i, int j, int k) const;
  double get_h() const { return _h; }

  virtual void update_impl(node * n, node ** nb, int parent, double & T) = 0;
  
EIKONAL_PRIVATE:
  void init();

  virtual void visit_neighbors_impl(abstract_node * n);

  // We initialize all of these variables to an invalid state so that
  // we can assert that we're using them correctly later if we've
  // compiled in debug mode.
  node * _nodes {nullptr};
  double const * _s_cache {nullptr};
  double _h {-1};
  int _height {-1}, _width {-1}, _depth {-1};
};

#include "marcher_3d.impl.hpp"

#endif // __MARCHER_3D_HPP__
