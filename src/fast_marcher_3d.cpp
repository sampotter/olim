#include "fast_marcher_3d.hpp"

#include <cassert>

void fast_marcher_3d::add_boundary_node(int i, int j, int k, double value) {
  assert(in_bounds(i, j, k));
  assert(this->operator()(i, j, k).is_far()); // TODO: for now---worried about heap
  this->operator()(i, j, k) = node_3d::make_boundary_node(i, j, k, value);
  stage_neighbors(i, j, k);
}

void fast_marcher_3d::stage_neighbor(int i, int j, int k) {
  if (in_bounds(i, j, k) && this->operator()(i, j, k).is_far()) {
    this->operator()(i, j, k).set_trial();
    insert_into_heap(&this->operator()(i, j, k));
  }
}

double fast_marcher_3d::get_value(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return this->operator()(i, j, k).get_value();
}

void fast_marcher_3d::update_node_value(int i, int j, int k) {
  assert(in_bounds(i, j, k));
  double T = std::numeric_limits<double>::infinity();
  update_node_value_impl(i, j, k, T);
  node_3d * n = &this->operator()(i, j, k);
  assert(n->is_trial());
  if (T <= n->get_value()) {
    n->set_value(T);
    adjust_heap_entry(n);
  }
}

node_3d & fast_marcher_3d::operator()(int i, int j, int k) {
  assert(in_bounds(i, j, k));
  return _nodes[_height*(_width*k + j) + i];
}

node_3d const & fast_marcher_3d::operator()(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return _nodes[_height*(_width*k + j) + i];
}

bool fast_marcher_3d::in_bounds(int i, int j, int k) const {
  return (unsigned) i < (unsigned) _height &&
    (unsigned) j < (unsigned) _width && (unsigned) k < (unsigned) _depth;
}

bool fast_marcher_3d::is_valid(int i, int j, int k) const {
  return in_bounds(i, j, k) && this->operator()(i, j, k).is_valid();
}

double fast_marcher_3d::S(int i, int j, int k) {
  assert(in_bounds(i, j, k));
  int l = _height*(_width*k + j) + i;
  assert(l < static_cast<int>(_S_cache.size()));
  if (_S_cache[l] < 0) {
    _S_cache[l] = _S(_h*j - _x0, _h*i - _y0, _h*k - _z0);
  }
  return _S_cache[l];
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
