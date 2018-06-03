#ifndef __ABSTRACT_MARCHER_HPP__
#define __ABSTRACT_MARCHER_HPP__

#include "heap.hpp"
#include "abstract_node.hpp"

struct abstract_marcher {
  void run();
  void step();
  virtual ~abstract_marcher() {}
EIKONAL_PROTECTED:
  abstract_marcher(size_t initial_heap_size = 256ul);
  void visit_neighbors(abstract_node * n);
  abstract_node * get_next_node();
  void adjust_heap_entry(abstract_node * n);
  void insert_into_heap(abstract_node * n);

  heap<abstract_node> _heap;
EIKONAL_PRIVATE:
  virtual void visit_neighbors_impl(abstract_node * n) = 0;
};

#endif // __ABSTRACT_MARCHER_HPP__
