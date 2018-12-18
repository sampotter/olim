#include "abstract_marcher.hpp"

#include <cmath>

void abstract_marcher::run() {
  abstract_node * n {nullptr};
  while (!_heap.empty()) {
    n = get_next_node();
    n->set_valid();
    visit_neighbors(n);
  }
}

void abstract_marcher::step() {
  auto * n = get_next_node();
  n->set_valid();
  visit_neighbors(n);
}

abstract_marcher::abstract_marcher(size_t initial_heap_size):
  _heap {initial_heap_size}
{}

void abstract_marcher::visit_neighbors(abstract_node * n) {
  visit_neighbors_impl(n);
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
