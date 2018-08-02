#ifndef __MARCHER_IMPL_HPP_HPP__
#define __MARCHER_IMPL_HPP_HPP__

#include <cassert>
#include <cmath>
#if PRINT_UPDATES
#    include <cstdio>
#endif // PRINT_UPDATES

#include "common.macros.hpp"
#include "offsets.hpp"

#define __di(k) di<2>[k]
#define __dj(k) dj<2>[k]

static inline size_t get_initial_heap_size(int width, int height) {
  return static_cast<size_t>(std::max(8.0, std::log(width*height)));
}

template <class base, class node>
marcher<base, node>::marcher(int height, int width, double h,
                             no_speed_func_t const &):
  abstract_marcher {get_initial_heap_size(width, height)},
  _nodes {new node[width*height]},
  _s_cache {new double[width*height]},
  _h {h},
  _height {height},
  _width {width}
{
  init();
}

template <class base, class node>
marcher<base, node>::marcher(int height, int width, double h,
                             double const * s_cache):
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

template <class base, class node>
marcher<base, node>::marcher(int height, int width, double h,
                       std::function<double(double, double)> s,
                       double x0, double y0):
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

template <class base, class node>
marcher<base, node>::~marcher()
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
template <class base, class node>
void marcher<base, node>::add_boundary_node(int i, int j, double value) {
#if PRINT_UPDATES
  printf("add_boundary_node(i = %d, j = %d, value = %g)\n", i, j, value);
#endif // PRINT_UPDATES
  assert(in_bounds(i, j));
  assert(operator()(i, j).is_far());
  visit_neighbors(&(operator()(i, j) = {i, j, value}));
}

template <class base, class node>
void marcher<base, node>::add_boundary_node(double x, double y, double value) {
  auto const dist = [x, y] (int i, int j) -> double {
    return std::sqrt((i - y)*(i - y) + (j - x)*(j - x));
  };

  int i0 = floor(y), i1 = ceil(y);
  int j0 = floor(x), j1 = ceil(x);

  node nodes[4] = {
    {i0, j0, value + dist(i0, j0)},
    {i0, j1, value + dist(i0, j1)},
    {i1, j0, value + dist(i1, j0)},
    {i1, j1, value + dist(i1, j1)}
  };

  add_boundary_nodes(nodes, 4);
}

template <class base, class node>
void marcher<base, node>::add_boundary_nodes(node const * nodes, int num) {
#if PRINT_UPDATES
  printf("add_boundary_nodes(nodes = %p, num = %d)\n", nodes, num);
#endif // PRINT_UPDATES

  node const * n;
  int i, j;

  /**
   * First, add the sequence of nodes to the grid.
   */
  for (int k = 0; k < num; ++k) {
    n = &nodes[k];
    i = n->get_i();
    j = n->get_j();
    assert(in_bounds(i, j));
    assert(operator()(i, j).is_far());
    operator()(i, j) = {i, j, n->get_value()};
  }

  /**
   * Next, with the nodes added, update their neighbors---this is done
   * in order to avoid breaking the heap property of the underlying
   * min-heap. Batch adding of boundary nodes may also be more
   * efficient for a large number of boundary nodes? (a guess)
   */
  for (int k = 0; k < num; ++k) {
    n = &nodes[k];
    i = n->get_i();
    j = n->get_j();
#if PRINT_UPDATES
    printf("add_boundary_node(i = %d, j = %d, value = %g)\n",
           i, j, n->get_value());
#endif // PRINT_UPDATES
    visit_neighbors(&operator()(i, j));
  }
}

template <class base, class node>
void marcher<base, node>::set_node_parent(int i, int j, int i_par, int j_par) {
  assert(in_bounds(i, j));
  assert(in_bounds(i_par, j_par));
  operator()(i, j).set_parent(&operator()(i_par, j_par));
}

template <class base, class node>
double marcher<base, node>::get_value(int i, int j) const {
  assert(in_bounds(i, j));
  return operator()(i, j).get_value();
}

template <class base, class node>
node & marcher<base, node>::operator()(int i, int j) {
  assert(in_bounds(i, j));
  assert(_nodes != nullptr);
  return _nodes[_width*i + j];
}

template <class base, class node>
node const & marcher<base, node>::operator()(int i, int j) const {
  assert(in_bounds(i, j));
  assert(_nodes != nullptr);
  return _nodes[_width*i + j];
}

template <class base, class node>
bool marcher<base, node>::in_bounds(int i, int j) const {
  return (unsigned) i < (unsigned) _height && (unsigned) j < (unsigned) _width;
}

template <class base, class node>
double marcher<base, node>::get_speed(int i, int j) const {
  assert(in_bounds(i, j));
  assert(_s_cache != nullptr);
  return _s_cache[_width*i + j];
}

template <class base, class node>
bool marcher<base, node>::is_valid(int i, int j) const {
  return in_bounds(i, j) && operator()(i, j).is_valid();
}

template <class base, class node>
void marcher<base, node>::init() {
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

template <class base, class node>
void marcher<base, node>::visit_neighbors_impl(abstract_node * n) {
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
#if PRINT_UPDATES
  printf("marcher::visit_neighbors_impl(i = %d, j = %d)\n",
         i, j);
#endif

  // These are temporary indices used below, analogous to i and j,
  // respectively.
  int a, b;

  // Traverse the update neighborhood of n and set all far nodes to
  // trial and insert them into the heap.
  for (int k = 0; k < base::nneib; ++k) {
    a = i + __di(k), b = j + __dj(k);
    if (in_bounds(a, b) && operator()(a, b).is_far()) {
      operator()(a, b).set_trial();
      insert_into_heap(&operator()(a, b));
    }
  }

  // Find the valid neighbors in the "full" neighborhood of n
  // (i.e. the unit max norm ball).
  node * valid_nb[8], * child_nb[base::nneib];
  memset(valid_nb, 0x0, 8*sizeof(abstract_node *));
  for (int k = 0; k < 8; ++k) {
    a = i + __di(k), b = j + __dj(k);
    if (in_bounds(a, b) && operator()(a, b).is_valid()) {
      valid_nb[k] = &this->operator()(a, b);
    }
  }

  // Some explanation of the indices used below:
  // - l is a radial index circling (a, b)
  // - parent is the radial index of (i, j) expressed in the same
  //   index space as l
  int di_k, dj_k;
  auto const set_child_nb = [&] (int parent) {
    memset(child_nb, 0x0, base::nneib*sizeof(abstract_node *));
    child_nb[parent] = static_cast<node *>(n);
    for (int l = 0; l < base::nneib; ++l) {
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
        child_nb[l] = valid_nb[m];
      }
    }
  };

  // Update the node at (i, j). Before calling, child_nb needs to be
  // filled appropriately. Upon updating, this sets the value of n and
  // adjusts its position in the heap.
  auto const update = [&] (int i, int j) {
    double T = INF(double);
    node * update_node = &operator()(i, j);
    update_impl(update_node, child_nb, T);
    if (T < update_node->get_value()) {
      update_node->set_value(T);
      adjust_heap_entry(update_node);
    }
  };

  // Get the parent index of a radial index `k'.
  auto const get_parent = [] (int k) {
    if (base::nneib == 4) {
      return (k + 2) % 4;
    } else {
      if (k < 4) return (k + 2) % 4;
      else return ((k - 2) % 4) + 4;
    }
  };

  // This the main update loop. Each neighbor of n which isn't valid
  // is now trial. For each neighboring trial node, use `set_child_nb'
  // to grab its valid neighbors and use `update' to actually update
  // the node's value and update its position in the heap.
  for (int k = 0; k < base::nneib; ++k) {
    if (!valid_nb[k]) {
      di_k = __di(k), dj_k = __dj(k);
      a = i + di_k, b = j + dj_k;
      if (!in_bounds(a, b)) continue;
      int parent = get_parent(k);
      set_child_nb(parent);
      update(a, b);
    }
  }
}

#undef __di
#undef __dj

#endif // __MARCHER_IMPL_HPP_HPP__
