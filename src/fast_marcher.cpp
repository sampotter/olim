#include "fast_marcher.hpp"

#include <cassert>
#include <cmath>
#include <vector>

fast_marcher::fast_marcher(size_t height, size_t width, double h, speed_func F):
  _nodes {new node[height*width]},
  _heap {static_cast<size_t>(std::log(height*width))}, // whatever
  _h {h},
  _F {F},
  _height {height},
  _width {width}
{
  for (size_t i = 0; i < _height; ++i) {
    for (size_t j = 0; j < _width; ++j) {
      this->operator()(i, j).set_i(i);
      this->operator()(i, j).set_j(j);
    }
  }
}

fast_marcher::~fast_marcher() {
  delete[] _nodes;
}

node & fast_marcher::operator()(size_t i, size_t j) {
  return _nodes[_width*i + j];
}

node const & fast_marcher::operator()(size_t i, size_t j) const {
  return _nodes[_width*i + j];
}

void fast_marcher::add_boundary_node(size_t i, size_t j) {
  this->operator()(i, j) = node::make_boundary_node(i, j);
  update_neighbors(i, j);
}

void fast_marcher::update_neighbors(size_t i, size_t j) {
  update_neighbors_impl(i, j);
}

void fast_marcher::run() {
  node* n = nullptr;
  while (!_heap.empty()) {
    n = get_next_node();
    n->set_valid();
    update_neighbors(n->get_i(), n->get_j());
  }
}

double fast_marcher::get_value(size_t i, size_t j) const {
  return this->operator()(i, j).get_value();
}

void fast_marcher::update_node_value(size_t i, size_t j) {
  update_node_value_impl(i, j);
}

bool fast_marcher::valid_index(size_t i, size_t j) const {
  return i < _height && j < _width;
}

node* fast_marcher::get_next_node() {
  auto const elt = _heap.front();
  _heap.pop_front();
  return elt;
}

double fast_marcher::get_h() const {
  return _h;
}

double fast_marcher::F(double x, double y) const {
  return _F(x, y);
}

double fast_marcher::F(size_t i, size_t j) const {
  // Use -1 to indicate that the cache value has not been set
  static std::vector<double> cache(_width*_height, -1);
  size_t k = _width*i + j;
  if (cache[k] < 0) {
    cache[k] = _F(_h*i, _h*j);
  }
  return cache[k];
}

void fast_marcher::adjust_heap_entry(node* n) {
  _heap.adjust_entry(n);
}

void fast_marcher::insert_into_heap(node* n) {
  _heap.insert(n);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
