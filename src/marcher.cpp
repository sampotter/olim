#include "marcher.hpp"

#include <cassert>
#include <vector>

marcher::marcher(int height, int width, double h):
  abstract_marcher {h, width*height},
  _nodes {new node[width*height]},
  _height {height},
  _width {width}
{
  init();
}

marcher::marcher(int height, int width, double h, speed_func S,
                       double x0, double y0):
  abstract_marcher {h, width*height},
  _nodes {new node[width*height]},
  _S {S},
  _x0 {x0},
  _y0 {y0},
  _height {height},
  _width {width}
{
  init();
}

marcher::marcher(int height, int width, double h, double * S_cache):
  abstract_marcher {h, width*height, S_cache},
  _nodes {new node[width*height]},
  _height {height},
  _width {width}
{
  init();
}

void marcher::add_boundary_node(int i, int j, double value) {
  assert(in_bounds(i, j));
  assert(this->operator()(i, j).is_far()); // TODO: for now---worried about heap
  this->operator()(i, j) = node::make_boundary_node(i, j, value);
  stage_neighbors(&this->operator()(i, j));
}

double marcher::get_value(int i, int j) const {
  assert(in_bounds(i, j));
  return this->operator()(i, j).get_value();
}

node & marcher::operator()(int i, int j) {
  assert(in_bounds(i, j));
  return _nodes[_width*i + j];
}

node const & marcher::operator()(int i, int j) const {
  assert(in_bounds(i, j));
  return _nodes[_width*i + j];
}

void marcher::update_node_value(int i, int j) {
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

void marcher::stage_neighbor(int i, int j) {
  if (in_bounds(i, j) && this->operator()(i, j).is_far()) {
    this->operator()(i, j).set_trial();
    insert_into_heap(&this->operator()(i, j));
  }
}

bool marcher::in_bounds(int i, int j) const {
  return (unsigned) i < (unsigned) _height && (unsigned) j < (unsigned) _width;
}

bool marcher::is_valid(int i, int j) const {
  return in_bounds(i, j) && this->operator()(i, j).is_valid();
}

double marcher::S(int i, int j) {
  assert(in_bounds(i, j));
  int k = _width*i + j;
  assert(k < static_cast<int>(_S_cache.size()));
  if (_S_cache[k] < 0) {
    _S_cache[k] = _S(_h*j - _x0, _h*i - _y0);
  }
  return _S_cache[k];
}

void marcher::init() {
  for (int i = 0; i < _height; ++i) {
    for (int j = 0; j < _width; ++j) {
      this->operator()(i, j).set_i(i);
      this->operator()(i, j).set_j(j);
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
