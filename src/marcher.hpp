#ifndef __MARCHER_HPP__
#define __MARCHER_HPP__

#include "abstract_marcher.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

#include <vector>

struct marcher: public abstract_marcher
{
  marcher(int height, int width, double h = 1,
          speed_func S = default_speed_func,
          double x0 = 0.0, double y0 = 0.0);
  marcher(int height, int width, double h, double * S_values);

  virtual void add_boundary_node(int i, int j, double value = 0.0) = 0;
  virtual double get_value(int i, int j) const = 0;

protected:
  bool in_bounds(int i, int j) const;
  double S(int i, int j);
  inline int get_height() const { return _height; }
  inline int get_width() const { return _width; }
  
  virtual void get_valid_neighbors(int i, int j, abstract_node ** nb) = 0;
  virtual void update_node_value_impl(int i, int j, double & T) = 0;

private:
  speed_func _S {default_speed_func};
  double _x0 {0}, _y0 {0};
  int _height;
  int _width;
};

#endif // __MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
