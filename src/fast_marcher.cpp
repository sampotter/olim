#include "fast_marcher.hpp"

#include <cassert>
#include <cmath>
#include <vector>

static size_t get_num_nodes(size_t height, size_t width) {
  return (height + 2)*(width + 2);
}

static size_t initial_heap_size(double height, double width) {
  return static_cast<size_t>(std::max(8.0, std::log(height*width)));
}

fast_marcher::fast_marcher(size_t height, size_t width, double h):
  _nodes {new node[get_num_nodes(height, width)]},
  _heap {initial_heap_size(height + 2, width + 2)},
  _h {h},
  _S_cache(get_num_nodes(height, width), -1),
  _height {height},
  _width {width}
{
  init();
}

fast_marcher::fast_marcher(size_t height, size_t width, double h, speed_func S,
                           double x0, double y0):
  _nodes {new node[get_num_nodes(height, width)]},
  _heap {initial_heap_size(height, width)},
  _h {h},
  _S_cache(get_num_nodes(height, width), -1),
  _S {S},
  _x0 {x0},
  _y0 {y0},
  _height {height},
  _width {width}
{
  init();
}

fast_marcher::fast_marcher(size_t height, size_t width, double h,
                           double const * const S_values):
  _nodes {new node[get_num_nodes(height, width)]},
  _heap {initial_heap_size(height, width)},
  _h {h},
  _S_cache(S_values, S_values + get_num_nodes(height, width)),
  _height {height},
  _width {width}
{
  init();
}

void fast_marcher::init() {
  for (size_t i = 0; i < _height + 2; ++i) {
    for (size_t j = 0; j < _width + 2; ++j) {
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
  return _nodes[get_linear_index(i, j)];
}

node const & fast_marcher::operator()(size_t i, size_t j) const {
  assert(in_bounds(i, j));
  return _nodes[get_linear_index(i, j)];
}

void fast_marcher::add_boundary_node(size_t i, size_t j) {
  size_t i_ = i + 1, j_ = j + 1;
  assert(in_bounds(i_, j_));
  this->operator()(i_, j_) = node::make_boundary_node(i_, j_);
  stage_neighbors(i_, j_);
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
  size_t i_ = i + 1, j_ = j + 1;
  assert(in_bounds(i_, j_));
  return this->operator()(i_, j_).get_value();
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
  return i < _height + 2 && j < _width + 2;
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

double fast_marcher::S(size_t i, size_t j) {
  assert(in_bounds(i, j));
  size_t k = get_linear_index(i, j);
  assert(k < _S_cache.size());
  if (_S_cache[k] < 0) {
    _S_cache[k] = _S(_h*(j + 1) - _x0, _h*(i + 1) - _y0);
  }
  return _S_cache[k];
}

void fast_marcher::adjust_heap_entry(node* n) {
  _heap.adjust_entry(n);
}

void fast_marcher::insert_into_heap(node* n) {
  _heap.insert(n);
}

size_t fast_marcher::get_linear_index(size_t i, size_t j) const {
  return (_width + 2)*i + j;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
