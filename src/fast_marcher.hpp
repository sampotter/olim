#ifndef __FAST_MARCHER_HPP__
#define __FAST_MARCHER_HPP__

#include "heap.hpp"
#include "node.hpp"
#include "typedefs.hpp"

#include <vector>

double default_speed_func(double x, double y) {
  (void) x;
  (void) y;
  return 1.0;
}

struct fast_marcher
{
  fast_marcher(size_t height, size_t width, double h = 1.0);
  fast_marcher(size_t height, size_t width, double h, speed_func S,
               double x0 = 0.0, double y0 = 0.0);
  fast_marcher(size_t height, size_t width, double h,
               double const * const S_values);
  virtual ~fast_marcher();

  void add_boundary_node(size_t i, size_t j);
  void run();
  double get_value(size_t i, size_t j) const;
    
protected:
  node & operator()(size_t i, size_t j);
  node const & operator()(size_t i, size_t j) const;
  void init();
  void update_node_value(size_t i, size_t j);
  void stage_neighbors(size_t i, size_t j);
  void stage_neighbor(size_t i, size_t j);
  virtual void get_valid_neighbors(size_t i, size_t j, node ** nb) = 0;
  bool in_bounds(size_t i, size_t j) const;
  bool is_valid(size_t i, size_t j) const;
  node* get_next_node();
  double get_h() const;
  double S(size_t i, size_t j);
  void adjust_heap_entry(node * n);
  void insert_into_heap(node * n);
	
private:
  virtual void update_node_value_impl(size_t i, size_t j, double & T) = 0;
  virtual void stage_neighbors_impl(size_t i, size_t j) = 0;

  node* _nodes;
  heap _heap;
  double _h;
  std::vector<double> _S_cache;
  speed_func _S {default_speed_func};
  double _x0 {0}, _y0 {0};
  size_t _height;
  size_t _width;
};

#endif // __FAST_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
