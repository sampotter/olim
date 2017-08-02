#ifndef __MARCHER_HPP__
#define __MARCHER_HPP__

#include "abstract_marcher.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

template <class Node>
struct marcher: public abstract_marcher {
  marcher(int height, int width, double h = 1,
          speed_func S = default_speed_func,
          double x0 = 0.0, double y0 = 0.0);
  marcher(int height, int width, double h, double * S_values);

  void add_boundary_node(int i, int j, double value = 0.0);
  double get_value(int i, int j) const;

protected:
  Node & operator()(int i, int j);
  Node const & operator()(int i, int j) const;
  void stage_neighbor(int i, int j);
  void update(int i, int j);
  bool in_bounds(int i, int j) const;
  double S(int i, int j);
  bool is_valid(int i, int j) const;
  
  virtual void get_valid_neighbors(int i, int j, abstract_node ** nb) = 0;
  virtual void update_impl(int i, int j, double & T) = 0;

private:
  void init();

  Node * _nodes;
  speed_func _S {default_speed_func};
  double _x0 {0}, _y0 {0};
  int _height;
  int _width;
};

#include "marcher.impl.hpp"

#endif // __MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
