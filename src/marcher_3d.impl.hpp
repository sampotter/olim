#ifndef __MARCHER_3D_IMPL_HPP__
#define __MARCHER_3D_IMPL_HPP__

#include <cassert>
#include <cmath>

#include "offsets.hpp"

#define __di(l) di<3>[l]
#define __dj(l) dj<3>[l]
#define __dk(l) dk<3>[l]

#define __linear_index(i, j, k) (_height*(_width*k + j) + i) // column-major
#define __x(l) (h*l - x0)
#define __y(l) (h*l - y0)
#define __z(l) (h*l - z0)

static inline size_t get_initial_heap_size(int width, int height, int depth) {
  return static_cast<size_t>(std::max(8.0, std::log(width*height*depth)));
}

template <class base, class node>
marcher_3d<base, node>::marcher_3d(int height, int width, int depth, double h,
                             no_speed_func_t const &):
  abstract_marcher {get_initial_heap_size(width, height, depth)},
  _nodes {new node[width*height*depth]},
  _s_cache {new double[width*height*depth]},
  _h {h},
  _height {height},
  _width {width},
  _depth {depth}
{
  init();
}

template <class base, class node>
marcher_3d<base, node>::marcher_3d(
  int height, int width, int depth, double h,
  std::function<double(double, double, double)> s,
  double x0, double y0, double z0):
  abstract_marcher {get_initial_heap_size(width, height, depth)},
  _nodes {new node[width*height*depth]},
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

template <class base, class node>
marcher_3d<base, node>::marcher_3d(int height, int width, int depth, double h,
                             double const * s_cache):
  abstract_marcher {get_initial_heap_size(width, height, depth)},
  _nodes {new node[width*height*depth]},
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

template <class base, class node>
marcher_3d<base, node>::~marcher_3d()
{
  delete[] _nodes;
  delete[] _s_cache;
}

/**
 * TODO: see comment about this function in marcher.impl.hpp (i.e. the
 * 2D version of this function).
 */
template <class base, class node>
void marcher_3d<base, node>::add_boundary_node(int i, int j, int k, double value) {
  assert(in_bounds(i, j, k));
  assert(operator()(i, j, k).is_far()); // TODO: for now---worried about heap
  visit_neighbors(&(operator()(i, j, k) = {i, j, k, value}));
}

template <class base, class node>
void marcher_3d<base, node>::add_boundary_node(double x, double y, double z,
                                         double value) {
  auto const dist = [x, y, z] (int i, int j, int k) -> double {
    return std::sqrt((i - y)*(i - y) + (j - x)*(j - x) + (k - z)*(k - z));
  };

  int i0 = floor(y), i1 = ceil(y);
  int j0 = floor(x), j1 = ceil(x);
  int k0 = floor(z), k1 = ceil(z);

  node nodes[8] = {
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

template <class base, class node>
void marcher_3d<base, node>::add_boundary_nodes(node const * nodes, int num) {
  node const * n;
  int i, j, k;

  /**
   * Add nodes to min-heap.
   */
  for (int l = 0; l < num; ++l) {
    n = &nodes[l];
    i = n->get_i();
    j = n->get_j();
    k = n->get_k();
    assert(in_bounds(i, j, k));
    assert(operator()(i, j, k).is_far());
    operator()(i, j, k) = {i, j, k, n->get_value()};
  }

  /**
   * Stage nodes' neighbors.
   */
  for (int l = 0; l < num; ++l) {
    n = &nodes[k];
    i = n->get_i();
    j = n->get_j();
    k = n->get_k();
    visit_neighbors(&operator()(i, j, k));
  }
}

template <class base, class node>
double marcher_3d<base, node>::get_value(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return operator()(i, j, k).get_value();
}

template <class base, class node>
node & marcher_3d<base, node>::operator()(int i, int j, int k) {
  assert(in_bounds(i, j, k));
  return _nodes[__linear_index(i, j, k)];
}

template <class base, class node>
node const & marcher_3d<base, node>::operator()(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return _nodes[__linear_index(i, j, k)];
}

template <class base, class node>
void marcher_3d<base, node>::stage(int i, int j, int k) {
  if (in_bounds(i, j, k) && operator()(i, j, k).is_far()) {
    operator()(i, j, k).set_trial();
    insert_into_heap(&operator()(i, j, k));
  }
}

template <class base, class node>
bool marcher_3d<base, node>::in_bounds(int i, int j, int k) const {
  return (unsigned) i < (unsigned) _height &&
    (unsigned) j < (unsigned) _width && (unsigned) k < (unsigned) _depth;
}

template <class base, class node>
double marcher_3d<base, node>::get_speed(int i, int j, int k) const {
  assert(in_bounds(i, j, k));
  return _s_cache[__linear_index(i, j, k)];
}

template <class base, class node>
bool marcher_3d<base, node>::is_valid(int i, int j, int k) const {
  return in_bounds(i, j, k) && operator()(i, j, k).is_valid();
}

template <class base, class node>
void marcher_3d<base, node>::init() {
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

template <class base, class node>
void marcher_3d<base, node>::visit_neighbors_impl(abstract_node * n) {
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
  int k = static_cast<node *>(n)->get_k();

#if PRINT_UPDATES
  printf("olim3d::visit_neighbors_impl(i = %d, j = %d, k = %d)\n",
         i, j, k);
#endif

  for (int l = 0; l < base::nneib; ++l) {
    this->stage(i + __di(l), j + __dj(l), k + __dk(l));
  }

  auto const update = [&] (int i, int j, int k, int parent) {
    double T = std::numeric_limits<double>::infinity();
    update_impl(i, j, k, parent, T);
    auto * n = &operator()(i, j, k);
    assert(n->is_trial());
    if (T < n->get_value()) {
      n->set_value(T);
      adjust_heap_entry(n);
    }
  };

  int a, b, c, l;
  for (l = 0; l < 6; ++l) {
    a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      update(a, b, c, (l + 3) % 6);
    }
  }
  if (base::nneib >= 18) {
    for (l = 6; l < 18; ++l) {
      a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
      if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
        update(a, b, c, 22 - 2*(l/2) + (l % 2));
      }
    }
  }
  if (base::nneib >= 26) {
    for (l = 18; l < 26; ++l) {
      a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
      if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
        update(a, b, c, 42 - 2*(l/2) + (l % 2));
      }
    }
  }
}

#undef __di
#undef __dj
#undef __dk

#undef __linear_index
#undef __x
#undef __y
#undef __z

#endif // __MARCHER_3D_IMPL_HPP__
