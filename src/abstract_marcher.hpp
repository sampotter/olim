#ifndef __ABSTRACT_MARCHER_HPP__
#define __ABSTRACT_MARCHER_HPP__

#include "heap.hpp"
#include "abstract_node.hpp"

struct abstract_marcher {
  void run();
protected:
  abstract_marcher(double h, int size, double * S_cache = nullptr);
  virtual ~abstract_marcher();
  void stage_neighbors(abstract_node * n);
  abstract_node * get_next_node();
  void adjust_heap_entry(abstract_node * n);
  void insert_into_heap(abstract_node * n);
  double get_h() const { return _h; }

  heap<abstract_node> _heap;
  double _h {1};
  double * _S_cache {nullptr};
private:
  virtual void stage_neighbors_impl(abstract_node * n) = 0;

  bool _should_free_S_cache;
};

#endif // __ABSTRACT_MARCHER_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
