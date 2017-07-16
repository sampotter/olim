#ifndef __MARCHER_2D_HPP__
#define __MARCHER_2D_HPP__

#include "abstract_marcher.hpp"
#include "node.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

#include <vector>

struct marcher_2d: public abstract_marcher
{
  marcher_2d(int height, int width, double h = 1.0);
  marcher_2d(int height, int width, double h, speed_func S,
             double x0 = 0.0, double y0 = 0.0);
  marcher_2d(int height, int width, double h, double * S_values);

  void add_boundary_node(int i, int j, double value = 0.0);
  double get_value(int i, int j) const;
    
protected:
  node & operator()(int i, int j);
  node const & operator()(int i, int j) const;
  void update_node_value(int i, int j);
  void stage_neighbor(int i, int j);
  bool in_bounds(int i, int j) const;
  bool is_valid(int i, int j) const;
  double S(int i, int j);
	
  virtual void get_valid_neighbors(int i, int j, node ** nb) = 0;
  
private:
  void init();
  
  virtual void update_node_value_impl(int i, int j, double & T) = 0;

  node * _nodes;
  speed_func _S {default_speed_func};
  double _x0 {0}, _y0 {0};
  int _height;
  int _width;
};

#endif // __MARCHER_2D_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
