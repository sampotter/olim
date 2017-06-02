#ifndef __HEAP_HPP__
#define __HEAP_HPP__

#include <cstddef>

#include "node.hpp"

struct heap
{
  heap(size_t capacity);
  ~heap();

  node* & front();
  node* const & front() const;
  bool empty() const;
  node** data() const;
  size_t size() const;
  void pop_front();
  void insert(node* n);
  void adjust_entry(node* n);
  void print() const;
private:
  void grow();
  void heapify(size_t pos);
  bool has_heap_property() const;
  void swap(size_t pos1, size_t pos2);
  double get_value(size_t pos) const;
	
  node** _data {nullptr};
  size_t _size {0};
  size_t _capacity {0};
};

#endif // __HEAP_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
