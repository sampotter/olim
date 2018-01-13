#ifndef __MARCHER_HPP__
#define __MARCHER_HPP__

#include "abstract_marcher.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

template <class Node>
struct marcher: public abstract_marcher {
  using float_type = double;
  using node_type = Node;

  static constexpr int ndim = 2;

  marcher(int height, int width, double h, no_speed_func_t const &);
  marcher(int height, int width, double h, double const * s_cache);
  marcher(int height, int width, double h = 1,
          std::function<double(double, double)> speed =
            static_cast<speed_func>(default_speed_func),
          double x0 = 0.0, double y0 = 0.0);
  virtual ~marcher();

  void add_boundary_node(int i, int j, double value = 0.0);
  void add_boundary_node(double i, double j, double value = 0.0);
  void add_boundary_nodes(Node const * nodes, int num_nodes);
  Node * get_node_pointer() const { return _nodes; }
  double get_value(int i, int j) const;
  int get_height() const { return _height; }
  int get_width() const { return _width; }
  void * get_s_cache_data() { return (void *) _s_cache; }

protected:
  Node & operator()(int i, int j);
  Node const & operator()(int i, int j) const;
  void stage(int i, int j);
  void update(int i, int j);
  bool in_bounds(int i, int j) const;
  double speed(int i, int j);
  bool is_valid(int i, int j) const;
  double get_h() const { return _h; }
  
  virtual void get_valid_neighbors(int i, int j, abstract_node ** nb) = 0;
  virtual void update_impl(int i, int j, double & T) = 0;

private:
  void init();

  Node * _nodes;
  double const * _s_cache {nullptr};
  double _h {1};
  int _height;
  int _width;
};

#include "marcher.impl.hpp"

#endif // __MARCHER_HPP__
