#ifndef __MARCHER_3D_IMPL_HPP__
#define __MARCHER_3D_IMPL_HPP__

#include <cassert>
#include <cmath>

#define __linear_index(i, j, k) (_height*(_width*k + j) + i) // column-major
#define __x(l) (h*l - x0)
#define __y(l) (h*l - y0)
#define __z(l) (h*l - z0)

static inline size_t get_initial_heap_size(int width, int height, int depth) {
  return static_cast<size_t>(std::max(8.0, std::log(width*height*depth)));
}

template <class Node>
marcher_3d<Node>::marcher_3d(int height, int width, int depth, double h,
                             no_speed_func_t const &):
  abstract_marcher {get_initial_heap_size(width, height, depth)},
  _nodes {new Node[width*height*depth]},
  _s_cache {new double[width*height*depth]},
  _h {h},
  _height {height},
  _width {width},
  _depth {depth}
{
  init();
}

template <class Node>
marcher_3d<Node>::marcher_3d(
  int height, int width, int depth, double h,
  std::function<double(double, double, double)> s,
  double x0, double y0, double z0):
  abstract_marcher {get_initial_heap_size(width, height, depth)},
  _nodes {new Node[width*height*depth]},
  _s_cache {new double[width*height*depth]},
  _h {h},
  _height {height},
  _width {width},
  _depth {depth}
{
  /**
   * Temporarily grab a writable pointer to cache the speed function
   * values.
   */
  // TODO: make sure this is being done in the most cache-friendly way
  // possible
  double * ptr = const_cast<double *>(_s_cache);
  for (int k = 0; k < depth; ++k) {
    for (int j = 0; j < width; ++j) {
      for (int i = 0; i < height; ++i) {
        ptr[__linear_index(i, j, k)] = s(__x(j), __y(i), __z(k));
      }
    }
  }

  init();
}

template <class Node>
marcher_3d<Node>::marcher_3d(int height, int width, int depth, double h,
                             double const * s_cache):
  abstract_marcher {get_initial_heap_size(width, height, depth)},
  _nodes {new Node[width*height*depth]},
  _s_cache {new double[width*height*depth]},
  _h {h},
  _height {height},
  _width {width},
  _depth {depth}
{
  memcpy((void *) _s_cache, (void *) s_cache,
         sizeof(double)*height*width*depth);
  init();
}

template <class Node>
marcher_3d<Node>::~marcher_3d()
{
  delete[] _nodes;
  delete[] _s_cache;
}

/**
 * TODO: see comment about this function in marcher.impl.hpp (i.e. the
 * 2D version of this function).
 */
template <class Node>
void marcher_3d<Node>::add_boundary_node(int i, int j, int k, double value) {
  assert(in_bounds(i, j, k));
  assert(operator()(i, j, k).is_far()); // TODO: for now---worried about heap
  stage_neighbors(&(operator()(i, j, k) = {i, j, k, value}));
}

template <class Node>
void marcher_3d<Node>::add_boundary_node(double x, double y, double z,
                                         double value) {
  auto const dist = [x, y, z] (int i, int j, int k) -> double {
    return std::sqrt((i - y)*(i - y) + (j - x)*(j - x) + (k - z)*(k - z));
  };

  int i0 = floor(y), i1 = ceil(y);
  int j0 = floor(x), j1 = ceil(x);
  int k0 = floor(z), k1 = ceil(z);

  Node nodes[8] = {
    {i0, j0, k0, value + dist(i0, j0, k0)},
    {i0, j0, k1, value + dist(i0, j0, k1)},
    {i0, j1, k0, value + dist(i0, j1, k0)},
    {i1, j0, k0, value + dist(i1, j0, k0)},
    {i0, j1, k1, value + dist(i0, j1, k1)},
    {i1, j0, k1, value + dist(i1, j0, k1)},
    {i1, j1, k0, value + dist(i1, j1, k0)},
    {i1, j1, k1, value + dist(i1, j1, k1)}
  };

  add_boundary_nodes(nodes, 8);
}

template <class Node>
void marcher_3d<Node>::add_boundary_nodes(Node const * nodes, int n) {
  Node const * node;
  int i, j, k;

  /**
   * Add nodes to min-heap.
   */
  for (int l = 0; l < n; ++l) {
    node = &nodes[l];
    i = node->get_i();
    j = node->get_j();
    k = node->get_k();
    assert(in_bounds(i, j, k));
    assert(operator()(i, j, k).is_far());
    operator()(i, j, k) = {i, j, k, node->get_value()};
  }

  /**
   * Stage nodes' neighbors.
   */
  for (int l = 0; l < n; ++l) {
    node = &nodes[k];
    i = node->get_i();
    j = node->get_j();
    k = node->get_k();
    stage_neighbors(&operator()(i, j, k));
  }
}

template <class Node>
double marcher_3d<Node>::get_value(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return operator()(i, j, k).get_value();
}

template <class Node>
Node & marcher_3d<Node>::operator()(int i, int j, int k) {
  assert(in_bounds(i, j, k));
  return _nodes[__linear_index(i, j, k)];
}

template <class Node>
Node const & marcher_3d<Node>::operator()(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return _nodes[__linear_index(i, j, k)];
}

template <class Node>
void marcher_3d<Node>::update(int i, int j, int k, int parent) {
  assert(in_bounds(i, j, k));
  double T = std::numeric_limits<double>::infinity();
  update_impl(i, j, k, parent, T);
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
double marcher_3d<Node>::get_speed(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return _s_cache[__linear_index(i, j, k)];
}

template <class Node>
bool marcher_3d<Node>::is_valid(int i, int j, int k) const {
  return in_bounds(i, j, k) && operator()(i, j, k).is_valid();
}

template <class Node>
void marcher_3d<Node>::init() {
  for (int i = 0; i < _height; ++i) {
    for (int j = 0; j < _width; ++j) {
      for (int k = 0; k < _depth; ++k) {
        operator()(i, j, k).set_i(i);
        operator()(i, j, k).set_j(j);
        operator()(i, j, k).set_k(k);
      }
    }
  }
}

#undef __linear_index
#undef __x
#undef __y
#undef __z

#endif // __MARCHER_3D_IMPL_HPP__
