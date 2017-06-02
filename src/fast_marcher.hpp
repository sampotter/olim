#ifndef __FAST_MARCHER_HPP__
#define __FAST_MARCHER_HPP__

#include "heap.hpp"
#include "node.hpp"
#include "typedefs.hpp"

struct fast_marcher
{
  fast_marcher(size_t height, size_t width, double h, speed_func F);
  virtual ~fast_marcher();

  node & operator()(size_t i, size_t j);
  node const & operator()(size_t i, size_t j) const;

  void add_boundary_node(size_t i, size_t j);
  void run();
  double get_value(size_t i, size_t j) const;
    
protected:
  void update_node_value(size_t i, size_t j);
  void update_neighbors(size_t i, size_t j);
  bool valid_index(size_t i, size_t j) const;
  node* get_next_node();
  double get_h() const;
  double F(double x, double y) const;
  double F(size_t i, size_t j) const;
  void adjust_heap_entry(node* n);
  void insert_into_heap(node* n);
	
private:
  virtual void update_node_value_impl(size_t i, size_t j) = 0;
  virtual void update_neighbors_impl(size_t i, size_t j) = 0;
  virtual void get_valid_neighbors(size_t i, size_t j, node ** nb) = 0;

  node* _nodes;
  heap _heap;
  double _h;
  speed_func _F;
  size_t _height;
  size_t _width;
};

#endif // __FAST_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
