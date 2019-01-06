#ifndef __MARCHER_3D_IMPL_HPP__
#define __MARCHER_3D_IMPL_HPP__

#include <src/config.hpp>

#include <assert.h>
#include <math.h>

#include "offsets.hpp"

#define __di(l) di<3>[l]
#define __dj(l) dj<3>[l]
#define __dk(l) dk<3>[l]

#define __linear_index(i, j, k) (_height*(_width*k + j) + i) // column-major
#define __x(l) (h*l - x0)
#define __y(l) (h*l - y0)
#define __z(l) (h*l - z0)

static inline size_t get_initial_heap_size(int width, int height, int depth) {
  return static_cast<size_t>(fmax(8.0, log(width*height*depth)));
}

template <class base, class node, int num_neighbors>
marcher_3d<base, node, num_neighbors>::marcher_3d() {}

template <class base, class node, int num_neighbors>
marcher_3d<base, node, num_neighbors>::marcher_3d(int height, int width, int depth, double h,
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

template <class base, class node, int num_neighbors>
marcher_3d<base, node, num_neighbors>::marcher_3d(
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
  // Grab a writable pointer to cache the speed function values.
  double x, z, * ptr = const_cast<double *>(_s_cache);
  for (int k = 0; k < depth; ++k) {
    z = __z(k);
    for (int j = 0; j < width; ++j) {
      x = __x(j);
      for (int i = 0; i < height; ++i) {
        ptr[__linear_index(i, j, k)] = s(x, __y(i), z);
      }
    }
  }

  init();
}

template <class base, class node, int num_neighbors>
marcher_3d<base, node, num_neighbors>::marcher_3d(int height, int width, int depth, double h,
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

template <class base, class node, int num_neighbors>
marcher_3d<base, node, num_neighbors>::~marcher_3d()
{
  assert(_nodes != nullptr);
  delete[] _nodes;

  assert(_s_cache != nullptr);
  delete[] _s_cache;
}

/**
 * TODO: see comment about this function in marcher.impl.hpp (i.e. the
 * 2D version of this function).
 */
template <class base, class node, int num_neighbors>
void marcher_3d<base, node, num_neighbors>::add_boundary_node(int i, int j, int k, double value) {
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
  assert(operator()(i, j, k).is_far());
#endif
  visit_neighbors(&(operator()(i, j, k) = {i, j, k, value}));
}

template <class base, class node, int num_neighbors>
void marcher_3d<base, node, num_neighbors>::add_boundary_node(double x, double y, double z,
                                               double value) {
  auto const dist = [x, y, z] (int i, int j, int k) -> double {
    return sqrt((i - y)*(i - y) + (j - x)*(j - x) + (k - z)*(k - z));
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

template <class base, class node, int num_neighbors>
void marcher_3d<base, node, num_neighbors>::add_boundary_nodes(node const * nodes, int num) {
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
#if EIKONAL_DEBUG && !RELWITHDEBINFO
    assert(in_bounds(i, j, k));
    assert(operator()(i, j, k).is_far());
#endif
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

template <class base, class node, int num_neighbors>
void marcher_3d<base, node, num_neighbors>::set_node_fac_parent(
  int i, int j, int k, int i_par, int j_par, int k_par)
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
  assert(in_bounds(i_par, j_par, k_par));
#endif
  operator()(i, j, k).set_fac_parent(&operator()(i_par, j_par, k_par));
}

template <class base, class node, int num_neighbors>
double marcher_3d<base, node, num_neighbors>::get_value(int i, int j, int k) const {
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
#endif
  return operator()(i, j, k).get_value();
}

template <class base, class node, int num_neighbors>
node & marcher_3d<base, node, num_neighbors>::operator()(int i, int j, int k) {
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
  assert(_nodes != nullptr);
#endif
  return _nodes[__linear_index(i, j, k)];
}

template <class base, class node, int num_neighbors>
node const & marcher_3d<base, node, num_neighbors>::operator()(int i, int j, int k) const {
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
  assert(_nodes != nullptr);
#endif
  return _nodes[__linear_index(i, j, k)];
}

template <class base, class node, int num_neighbors>
bool marcher_3d<base, node, num_neighbors>::in_bounds(int i, int j, int k) const {
  return (unsigned) i < (unsigned) _height &&
    (unsigned) j < (unsigned) _width && (unsigned) k < (unsigned) _depth;
}

template <class base, class node, int num_neighbors>
double marcher_3d<base, node, num_neighbors>::get_speed(int i, int j, int k) const {
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
  assert(_s_cache != nullptr);
#endif
  return _s_cache[__linear_index(i, j, k)];
}

template <class base, class node, int num_neighbors>
bool marcher_3d<base, node, num_neighbors>::is_valid(int i, int j, int k) const {
  return in_bounds(i, j, k) && operator()(i, j, k).is_valid();
}

template <class base, class node, int num_neighbors>
void marcher_3d<base, node, num_neighbors>::init() {
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

#define __maxabs3(x, y, z) \
  std::max(std::abs(x), std::max(std::abs(y), std::abs(z)))

template <class base, class node, int num_neighbors>
void marcher_3d<base, node, num_neighbors>::visit_neighbors_impl(abstract_node * n) {
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
  int k = static_cast<node *>(n)->get_k();

  // See comments in marcher.impl.hpp; the visit_neighbors_impl there
  // is done analogously to this one.

  int a, b, c;

  // Stage neighbors.
  for (int l = 0; l < num_neighbors; ++l) {
    a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
    if (in_bounds(a, b, c) && operator()(a, b, c).is_far()) {
      operator()(a, b, c).set_trial();
      insert_into_heap(&operator()(a, b, c));
    }
  }

  // Get valid neighbors.
  node * valid_nb[26], * child_nb[num_neighbors];
  memset(valid_nb, 0x0, 26*sizeof(abstract_node *));
  for (int l = 0; l < 26; ++l) {
    a = i + __di(l), b = j + __dj(l), c = k + __dk(l);
    if (in_bounds(a, b, c) && operator()(a, b, c).is_valid()) {
      valid_nb[l] = &this->operator()(a, b, c);
    }
  }

  int di_l, dj_l, dk_l;
  auto const set_child_nb = [&] (int parent) {
    memset(child_nb, 0x0, num_neighbors*sizeof(abstract_node *));
    child_nb[parent] = static_cast<node *>(n);
    for (int m = 0; m < num_neighbors; ++m) {
      if (m == parent) {
        continue;
      }
      int di_lm, dj_lm, dk_lm;
      if (std::abs(di_lm = di_l + __di(m)) > 1 ||
          std::abs(dj_lm = dj_l + __dj(m)) > 1 ||
          std::abs(dk_lm = dk_l + __dk(m)) > 1) {
        continue;
      }
      if (in_bounds(i + di_lm, j + dj_lm, k + dk_lm)) {
        int p = d2l(di_lm, dj_lm, dk_lm);
        child_nb[m] = valid_nb[p];
      }
    }
  };

  auto & s_hat = static_cast<base *>(this)->s_hat;
  auto const update = [&] (int i, int j, int k, int parent) {
    auto T = inf<double>;
    node * update_node = &operator()(i, j, k);
    s_hat = this->get_speed(i, j, k);
    update_impl(update_node, child_nb, parent, T);
#if NODE_MONITORING
    if (update_node->monitoring_node()) {
      std::cout << *update_node << std::endl;
    }
#endif
    if (T < update_node->get_value()) {
      update_node->set_value(T);
      adjust_heap_entry(update_node);
    }
  };

  auto const get_parent = [] (int l) {
    // TODO: check base::num_nb to reduce amount of branching here
    if (l < 6) return (l + 3) % 6;
    else if (l < 18) return 22 - 2*(l/2) + (l % 2);
    else return 42 - 2*(l/2) + (l % 2);
  };

  for (int l = 0; l < num_neighbors; ++l) {
    if (!valid_nb[l]) {
      di_l = __di(l), dj_l = __dj(l), dk_l = __dk(l);
      a = i + di_l, b = j + dj_l, c = k + dk_l;
      if (!in_bounds(a, b, c)) continue;
      int parent = get_parent(l);
      set_child_nb(parent);
      update(a, b, c, parent);
    }
  }
}

#undef __maxabs3

#undef __di
#undef __dj
#undef __dk

#undef __linear_index
#undef __x
#undef __y
#undef __z

#endif // __MARCHER_3D_IMPL_HPP__
