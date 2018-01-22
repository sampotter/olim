#include "abstract_marcher.hpp"

#include <src/config.hpp>
#if PRINT_UPDATES
#    include <cstdio>
#endif // PRINT_UPDATES

#include <cmath>

void abstract_marcher::run() {
#if PRINT_UPDATES
  puts("abstract_marcher::run()");
#endif // PRINT_UPDATES
  abstract_node * n {nullptr};
  while (!_heap.empty()) {
    n = get_next_node();
    n->set_valid();
    stage_neighbors(n);
  }
}

abstract_marcher::abstract_marcher(size_t initial_heap_size):
  _heap {initial_heap_size}
{}

void abstract_marcher::stage_neighbors(abstract_node * n) {
  stage_neighbors_impl(n);
}

abstract_node * abstract_marcher::get_next_node() {
  auto const n = _heap.front();
  _heap.pop_front();
  return n;
}

void abstract_marcher::adjust_heap_entry(abstract_node * n) {
  _heap.swim(n);
}

void abstract_marcher::insert_into_heap(abstract_node * n) {
  _heap.insert(n);
}
