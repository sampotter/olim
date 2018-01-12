#ifndef __ABSTRACT_MARCHER_HPP__
#define __ABSTRACT_MARCHER_HPP__

#include "heap.hpp"
#include "abstract_node.hpp"

struct abstract_marcher {
  void run();
  virtual ~abstract_marcher() {}
protected:
  abstract_marcher(size_t initial_heap_size = 256ul);
  void stage_neighbors(abstract_node * n);
  abstract_node * get_next_node();
  void adjust_heap_entry(abstract_node * n);
  void insert_into_heap(abstract_node * n);

  heap<abstract_node> _heap;
private:
  virtual void stage_neighbors_impl(abstract_node * n) = 0;
};

#endif // __ABSTRACT_MARCHER_HPP__
