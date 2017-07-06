#include "fast_marcher.hpp"

#include <cassert>
#include <cmath>
#include <vector>

static size_t initial_heap_size(int height, int width) {
  return static_cast<size_t>(std::max(8.0, std::log(height*width)));
}

fast_marcher::fast_marcher(int height, int width, double h):
  _nodes {new node[height*width]},
  _heap {initial_heap_size(height, width)},
  _h {h},
  _S_cache(height*width, -1),
  _height {height},
  _width {width}
{
  init();
}

fast_marcher::fast_marcher(int height, int width, double h, speed_func S,
                           double x0, double y0):
  _nodes {new node[height*width]},
  _heap {initial_heap_size(height, width)},
  _h {h},
  _S_cache(height*width, -1),
  _S {S},
  _x0 {x0},
  _y0 {y0},
  _height {height},
  _width {width}
{
  init();
}

fast_marcher::fast_marcher(int height, int width, double h,
                           double const * const S_values):
  _nodes {new node[height*width]},
  _heap {initial_heap_size(height, width)},
  _h {h},
  _S_cache(S_values, S_values + height*width),
  _height {height},
  _width {width}
{
  init();
}

void fast_marcher::init() {
  for (int i = 0; i < _height; ++i) {
    for (int j = 0; j < _width; ++j) {
      this->operator()(i, j).set_i(i);
      this->operator()(i, j).set_j(j);
    }
  }
}

fast_marcher::~fast_marcher() {
  delete[] _nodes;
}

node & fast_marcher::operator()(int i, int j) {
  assert(in_bounds(i, j));
  return _nodes[get_linear_index(i, j)];
}

node const & fast_marcher::operator()(int i, int j) const {
  assert(in_bounds(i, j));
  return _nodes[get_linear_index(i, j)];
}

void fast_marcher::add_boundary_node(int i, int j, double value) {
  assert(in_bounds(i, j));
  assert(this->operator()(i, j).is_far()); // TODO: for now---worried about heap
  this->operator()(i, j) = node::make_boundary_node(i, j, value);
  stage_neighbors(i, j);
}

void fast_marcher::stage_neighbors(int i, int j) {
  assert(in_bounds(i, j));
  stage_neighbors_impl(i, j);
}

void fast_marcher::stage_neighbor(int i, int j) {
  if (in_bounds(i, j) && this->operator()(i, j).is_far()) {
    this->operator()(i, j).set_trial();
    insert_into_heap(&this->operator()(i, j));
  }
}

void fast_marcher::run() {
  node * n = nullptr;
  while (!_heap.empty()) {
    n = get_next_node();
    n->set_valid();
    stage_neighbors(n->get_i(), n->get_j());
  }
}

double fast_marcher::get_value(int i, int j) const {
  assert(in_bounds(i, j));
  return this->operator()(i, j).get_value();
}

void fast_marcher::update_node_value(int i, int j) {
  assert(in_bounds(i, j));
  double T = std::numeric_limits<double>::infinity();
  update_node_value_impl(i, j, T);
  node * n = &this->operator()(i, j);
  assert(n->is_trial());
  if (T <= n->get_value()) {
    n->set_value(T);
    adjust_heap_entry(n);
  }
}

bool fast_marcher::in_bounds(int i, int j) const {
  return 0 <= i && i < _height && 0 <= j && j < _width;
}

bool fast_marcher::is_valid(int i, int j) const {
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

double fast_marcher::S(int i, int j) {
  assert(in_bounds(i, j));
  int k = get_linear_index(i, j);
  assert(k < static_cast<int>(_S_cache.size()));
  if (_S_cache[k] < 0) {
    _S_cache[k] = _S(_h*j - _x0, _h*i - _y0);
  }
  return _S_cache[k];
}

void fast_marcher::adjust_heap_entry(node* n) {
  _heap.adjust_entry(n);
}

void fast_marcher::insert_into_heap(node* n) {
  _heap.insert(n);
}

int fast_marcher::get_linear_index(int i, int j) const {
  return _width*i + j;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
