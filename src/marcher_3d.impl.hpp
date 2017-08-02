#ifndef __MARCHER_3D_IMPL_HPP__
#define __MARCHER_3D_IMPL_HPP__

#include <cassert>

template <class Node>
void marcher_3d<Node>::add_boundary_node(int i, int j, int k, double value) {
  assert(in_bounds(i, j, k));
  assert(operator()(i, j, k).is_far()); // TODO: for now---worried about heap
  stage_neighbors(&(operator()(i, j, k) = {i, j, k, value}));
}

template <class Node>
double marcher_3d<Node>::get_value(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return operator()(i, j, k).get_value();
}

template <class Node>
Node & marcher_3d<Node>::operator()(int i, int j, int k) {
  assert(in_bounds(i, j, k));
  return _nodes[_height*(_width*k + j) + i];
}

template <class Node>
Node const & marcher_3d<Node>::operator()(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return _nodes[_height*(_width*k + j) + i];
}

template <class Node>
void marcher_3d<Node>::update(int i, int j, int k) {
  assert(in_bounds(i, j, k));
  double T = std::numeric_limits<double>::infinity();
  update_impl(i, j, k, T);
  auto * n = &operator()(i, j, k);
  assert(n->is_trial());
  if (T <= n->get_value()) {
    n->set_value(T);
    adjust_heap_entry(n);
  }
}

template <class Node>
void marcher_3d<Node>::stage(int i, int j, int k) {
  if (in_bounds(i, j, k) && operator()(i, j, k).is_far()) {
    operator()(i, j, k).set_trial();
    insert_into_heap(&operator()(i, j, k));
  }
}

template <class Node>
bool marcher_3d<Node>::in_bounds(int i, int j, int k) const {
  return (unsigned) i < (unsigned) _height &&
    (unsigned) j < (unsigned) _width && (unsigned) k < (unsigned) _depth;
}

template <class Node>
bool marcher_3d<Node>::is_valid(int i, int j, int k) const {
  return in_bounds(i, j, k) && operator()(i, j, k).is_valid();
}

template <class Node>
double marcher_3d<Node>::S(int i, int j, int k) {
  assert(in_bounds(i, j, k));
  int l = _height*(_width*k + j) + i;
  assert(l < static_cast<int>(_S_cache.size()));
  if (_S_cache[l] < 0) {
    _S_cache[l] = _S(_h*j - _x0, _h*i - _y0, _h*k - _z0);
  }
  return _S_cache[l];
}

#endif // __MARCHER_3D_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
