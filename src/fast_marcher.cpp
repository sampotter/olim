#include "fast_marcher.hpp"

#include <cassert>
#include <cmath>
#include <vector>

static double default_speed_func(double x, double y) {
  (void) x;
  (void) y;
  return 1.0;
}

fast_marcher::fast_marcher(size_t height, size_t width, double h):
  _nodes {new node[height*width]},
  _heap {static_cast<size_t>(std::max(8.0, std::log(height*width)))}, // whatever
  _h {h},
  _F_cache(width*height, -1),
  _F {default_speed_func},
  _height {height},
  _width {width}
{
  init();
}

fast_marcher::fast_marcher(size_t height, size_t width, double h, speed_func F):
  _nodes {new node[height*width]},
  _heap {static_cast<size_t>(std::log(height*width))}, // whatever
  _h {h},
  _F_cache(width*height, -1),
  _F {F},
  _height {height},
  _width {width}
{
  init();
}

void fast_marcher::init() {
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
  assert(in_bounds(i, j));
  return _nodes[_width*i + j];
}

node const & fast_marcher::operator()(size_t i, size_t j) const {
  assert(in_bounds(i, j));
  return _nodes[_width*i + j];
}

void fast_marcher::add_boundary_node(size_t i, size_t j) {
  assert(in_bounds(i, j));
  this->operator()(i, j) = node::make_boundary_node(i, j);
  stage_neighbors(i, j);
}

void fast_marcher::stage_neighbors(size_t i, size_t j) {
  assert(in_bounds(i, j));
  stage_neighbors_impl(i, j);
}

void fast_marcher::stage_neighbor(size_t i, size_t j) {
  if (in_bounds(i, j) && this->operator()(i, j).is_far()) {
    this->operator()(i, j).set_trial();
    insert_into_heap(&this->operator()(i, j));
  }
}

void fast_marcher::run() {
  node* n = nullptr;
  while (!_heap.empty()) {
    n = get_next_node();
    n->set_valid();
    stage_neighbors(n->get_i(), n->get_j());
  }
}

double fast_marcher::get_value(size_t i, size_t j) const {
  assert(in_bounds(i, j));
  return this->operator()(i, j).get_value();
}

void fast_marcher::update_node_value(size_t i, size_t j) {
  assert(in_bounds(i, j));
  double T = std::numeric_limits<double>::infinity();
  update_node_value_impl(i, j, T);
  node* n = &this->operator()(i, j);
  assert(n->is_trial());
  if (T <= n->get_value()) {
    n->set_value(T);
    adjust_heap_entry(n);
  }
}

bool fast_marcher::in_bounds(size_t i, size_t j) const {
  return i < _height && j < _width;
}

bool fast_marcher::is_valid(size_t i, size_t j) const {
  return in_bounds(i, j) && this->operator()(i, j).is_valid();
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

double fast_marcher::F(size_t i, size_t j) {
  // Use -1 to indicate that the cache value has not been set
  assert(in_bounds(i, j));
  size_t k = _width*i + j;
  assert(k < _F_cache.size());
  if (_F_cache[k] < 0) {
    _F_cache[k] = _F(_h*i, _h*j);
  }
  return _F_cache[k];
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
