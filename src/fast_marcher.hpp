#ifndef __FAST_MARCHER_HPP__
#define __FAST_MARCHER_HPP__

#include "abstract_marcher.hpp"
#include "heap.hpp"
#include "node.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

#include <vector>

struct fast_marcher: public abstract_marcher
{
  fast_marcher(int height, int width, double h = 1.0);
  fast_marcher(int height, int width, double h, speed_func S,
               double x0 = 0.0, double y0 = 0.0);
  fast_marcher(int height, int width, double h, double * S_values);

  void add_boundary_node(int i, int j, double value = 0.0);
  void run();
  double get_value(int i, int j) const;
    
protected:
  node & operator()(int i, int j);
  node const & operator()(int i, int j) const;
  void init();
  void update_node_value(int i, int j);
  void stage_neighbors(int i, int j);
  void stage_neighbor(int i, int j);
  virtual void get_valid_neighbors(int i, int j, node ** nb) = 0;
  bool in_bounds(int i, int j) const;
  bool is_valid(int i, int j) const;
  node * get_next_node();
  double S(int i, int j);
  void adjust_heap_entry(node * n);
  void insert_into_heap(node * n);
	
private:
  virtual void update_node_value_impl(int i, int j, double & T) = 0;
  virtual void stage_neighbors_impl(int i, int j) = 0;

  node * _nodes;
  heap<node> _heap;
  speed_func _S {default_speed_func};
  double _x0 {0}, _y0 {0};
  int _height;
  int _width;
};

#endif // __FAST_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
