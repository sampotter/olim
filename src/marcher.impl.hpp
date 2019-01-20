#ifndef __MARCHER_IMPL_HPP_HPP__
#define __MARCHER_IMPL_HPP_HPP__

#include <assert.h>
#include <math.h>

#include "common.hpp"
#include "offsets.hpp"
#include "updates.line.hpp"

#define __di(k) di<2>[k]
#define __dj(k) dj<2>[k]

static inline size_t get_initial_heap_size(int width, int height) {
  return static_cast<size_t>(fmax(8.0, log(width*height)));
}

template <class base, class node, int num_neighbors>
marcher<base, node, num_neighbors>::marcher(
  int height, int width, double h, no_speed_func_t const &):
  abstract_marcher {get_initial_heap_size(width, height)},
  _nodes {new node[width*height]},
  _s_cache {new double[width*height]},
  _h {h},
  _height {height},
  _width {width}
{
  init();
}

template <class base, class node, int num_neighbors>
marcher<base, node, num_neighbors>::marcher(
  int height, int width, double h, double const * s_cache):
  abstract_marcher {get_initial_heap_size(width, height)},
  _nodes {new node[width*height]},
  _s_cache {new double[width*height]},
  _h {h},
  _height {height},
  _width {width}
{
  memcpy((void *) _s_cache, (void *) s_cache, sizeof(double)*height*width);
  init();
}

template <class base, class node, int num_neighbors>
marcher<base, node, num_neighbors>::marcher(
  int height, int width, double h,
  std::function<double(double, double)> s, double x0, double y0):
  abstract_marcher {get_initial_heap_size(width, height)},
  _nodes {new node[width*height]},
  _s_cache {new double[width*height]},
  _h {h},
  _height {height},
  _width {width}
{
  double * ptr = const_cast<double *>(_s_cache);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      ptr[i*width + j] = s(h*j - x0, h*i - y0);
    }
  }
  init();
}

template <class base, class node, int num_neighbors>
marcher<base, node, num_neighbors>::~marcher()
{
  assert(_nodes != nullptr);
  delete[] _nodes;

  assert(_s_cache != nullptr);
  delete[] _s_cache;
}

/**
 * TODO: this function is a little broken right now. If we remove the
 * "is_far" assertion and allow the user to add boundary nodes
 * wherever, then it's possible to invalidate the heap property in the
 * underlying min-heap. One way to get this to happen is to add a
 * boundary node and then add another boundary node which is adjacent
 * to the node just added. This should be fixed, since there is a
 * correct way to enable this behavior.
 */
template <class base, class node, int num_neighbors>
void
marcher<base, node, num_neighbors>::add_boundary_node(int i, int j, double value)
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
  assert(operator()(i, j).is_far());
#endif
  visit_neighbors(&(operator()(i, j) = {i, j, value}));
}

#define LINE(p0, u0, s, s0, h)                          \
  updates::line<base::F_>()(norm2<2>(p0), u0, s, s0, h)

template <class base, class node, int num_neighbors>
void
marcher<base, node, num_neighbors>::add_boundary_node(
  double x, double y, double s, double value)
{
#if PRINT_UPDATES
  printf("add_boundary_node(x = %g, y = %g, s = %g, value = %g)\n",
         x, y, s, value);
#endif
  double h = get_h(), i = y/h, j = x/h, u0 = value, s0 = s;
  assert(in_bounds(i, j));

  // TODO: this isn't as general as it could be. We also want to
  // handle cases where i or j are grid-aligned, in which case we need
  // to process 6 nodes. If i and j are both grid aligned, then we
  // should just call the integer add_boundary_node.
  int is[2] = {(int) floor(i), (int) floor(i) + 1};
  int js[2] = {(int) floor(j), (int) floor(j) + 1};

  double P[4][2] = {
    {i - is[0], j - js[0]},
    {i - is[1], j - js[0]},
    {i - is[0], j - js[1]},
    {i - is[1], j - js[1]}
  };

  for (int a = 0; a < 4; ++a) {
    int b0 = a & 1, b1 = (a & 2) >> 1;
    int i_ = is[b0], j_ = js[b1];
    assert(in_bounds(i_, j_));
    assert(operator()(i_, j_).is_far());
    double s_hat = get_speed(i_, j_), u_hat = LINE(P[a], u0, s_hat, s0, h);
    insert_into_heap(&(operator()(i_, j_) = {i_, j_, u_hat, state::trial}));
  }
}

#undef LINE

template <class base, class node, int num_neighbors>
void
marcher<base, node, num_neighbors>::add_boundary_nodes(
  node const * nodes, int num)
{
  node const * const * tmp = malloc(sizeof(node const * const *)*num);
  for (int i = 0; i < num; ++i) {
    tmp[i] = &nodes[i];
  }
  add_boundary_nodes(tmp, num);
  free(tmp);
}

template <class base, class node, int num_neighbors>
void
marcher<base, node, num_neighbors>::add_boundary_nodes(
  node const * const * nodes, int num)
{
  for (int k = 0; k < num; ++k) {
    auto n = nodes[k];
    int i = n->get_i(), j = n->get_j();
    assert(in_bounds(i, j));
    assert(operator()(i, j).is_far());
    double u = n->get_value();
    insert_into_heap(&(operator()(i, j) = {i, j, u, state::trial}));
  }
}

template <class base, class node, int num_neighbors>
void
marcher<base, node, num_neighbors>::set_node_fac_center(
  int i, int j, typename node::fac_center const * fc)
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
  assert(in_bounds(fc->i, fc->j));
#endif
  operator()(i, j).set_fac_center(fc);
}

template <class base, class node, int num_neighbors>
double
marcher<base, node, num_neighbors>::get_value(int i, int j) const
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
#endif
  return operator()(i, j).get_value();
}

template <class base, class node, int num_neighbors>
node &
marcher<base, node, num_neighbors>::operator()(int i, int j)
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
  assert(_nodes != nullptr);
#endif
  return _nodes[_width*i + j];
}

template <class base, class node, int num_neighbors>
node const &
marcher<base, node, num_neighbors>::operator()(int i, int j) const
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
  assert(_nodes != nullptr);
#endif
  return _nodes[_width*i + j];
}

template <class base, class node, int num_neighbors>
bool
marcher<base, node, num_neighbors>::in_bounds(int i, int j) const
{
  return (unsigned) i < (unsigned) _height && (unsigned) j < (unsigned) _width;
}

template <class base, class node, int num_neighbors>
bool marcher<base, node, num_neighbors>::in_bounds(double i, double j) const {
  return 0 <= i <= _height - 1 && 0 <= j <= _width - 1;
}

template <class base, class node, int num_neighbors>
double
marcher<base, node, num_neighbors>::get_speed(int i, int j) const
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
  assert(_s_cache != nullptr);
#endif
  return _s_cache[_width*i + j];
}

template <class base, class node, int num_neighbors>
bool
marcher<base, node, num_neighbors>::is_valid(int i, int j) const
{
  return in_bounds(i, j) && operator()(i, j).is_valid();
}

template <class base, class node, int num_neighbors>
void
marcher<base, node, num_neighbors>::init()
{
  /**
   * Set the indices associated with each node in the grid.
   */
  for (int i = 0; i < _height; ++i) {
    for (int j = 0; j < _width; ++j) {
      operator()(i, j).set_i(i);
      operator()(i, j).set_j(j);
    }
  }
}

template <class base, class node, int num_neighbors>
void
marcher<base, node, num_neighbors>::visit_neighbors_impl(abstract_node * n)
{
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();

  // These are temporary indices used below, analogous to i and j,
  // respectively.
  int a, b;

  // Traverse the update neighborhood of n and set all far nodes to
  // trial and insert them into the heap.
  for (int k = 0; k < num_neighbors; ++k) {
    a = i + __di(k), b = j + __dj(k);
    if (in_bounds(a, b) && operator()(a, b).is_far()) {
      operator()(a, b).set_trial();
      insert_into_heap(&operator()(a, b));
    }
  }

  // Find the valid neighbors in the "full" neighborhood of n
  // (i.e. the unit max norm ball).
  memset(valid, 0x0, 8*sizeof(abstract_node *));
  for (int k = 0; k < 8; ++k) {
    a = i + __di(k), b = j + __dj(k);
    if (in_bounds(a, b) && operator()(a, b).is_valid()) {
      valid[k] = &this->operator()(a, b);
    }
  }

  // Some explanation of the indices used below:
  // - l is a radial index circling (a, b)
  // - parent is the radial index of (i, j) expressed in the same
  //   index space as l
  int di_k, dj_k;
  node ** nb = static_cast<base *>(this)->nb;
  auto const set_nb = [&] (int parent) {
    memset(nb, 0x0, num_neighbors*sizeof(abstract_node *));
    nb[parent] = static_cast<node *>(n);
    for (int l = 0; l < num_neighbors; ++l) {
      if (l == parent) {
        continue;
      }
      int di_kl, dj_kl;
      if (std::abs(di_kl = di_k + __di(l)) > 1 ||
          std::abs(dj_kl = dj_k + __dj(l)) > 1) {
        continue;
      }
      if (in_bounds(i + di_kl, j + dj_kl)) {
        int m = d2l(di_kl, dj_kl);
        nb[l] = valid[m];
      }
    }
  };

  // Update the node at (i, j). Before calling, `nb' needs to be
  // filled appropriately. Upon updating, this sets the value of n and
  // adjusts its position in the heap.
  auto const update = [&] (int i, int j) {
    auto T = inf<double>;
    node * update_node = &operator()(i, j);
    static_cast<base *>(this)->s_hat = this->get_speed(i, j);
    update_impl(update_node, T);
    if (T < update_node->get_value()) {
      update_node->set_value(T);
      adjust_heap_entry(update_node);
    }
  };

  // Get the parent index of a radial index `k'.
  auto const get_parent = [] (int k) {
    if (num_neighbors == 4) {
      return (k + 2) % 4;
    } else {
      if (k < 4) return (k + 2) % 4;
      else return ((k - 2) % 4) + 4;
    }
  };

  // This is the main update loop. Each neighbor of n which isn't
  // valid is now trial. For each neighboring trial node, use
  // `set_child_nb' to grab its valid neighbors and use `update' to
  // actually update the node's value and update its position in the
  // heap.
  for (int k = 0; k < num_neighbors; ++k) {
    if (!valid[k]) {
      di_k = __di(k), dj_k = __dj(k);
      a = i + di_k, b = j + dj_k;
      if (!in_bounds(a, b)) continue;
      int parent = get_parent(k);
      set_nb(parent);
      update(a, b);
    }
  }
}

#undef __di
#undef __dj

#endif // __MARCHER_IMPL_HPP_HPP__
