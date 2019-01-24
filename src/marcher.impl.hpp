#ifndef __MARCHER_IMPL_HPP_HPP__
#define __MARCHER_IMPL_HPP_HPP__

#include <assert.h>
#include <math.h>

#include "common.hpp"
#include "offsets.hpp"
#include "updates.line.hpp"

#define __di(k) di<2>[k]
#define __dj(k) dj<2>[k]

// TODO: really need an external memory constructor that will let us
// use external memory somewhere for _s_cache and _U so we don't have
// to double up on these...

template <class base, int num_nb>
marcher<base, num_nb>::marcher(
  int height, int width, double h, no_speed_func_t const &):
  _size {width*height},
  _heap {{this}, initial_heap_capacity(_size)},
  _U {new double[_size]},
  _s_cache {new double[_size]},
  _state {new state[_size]},
  _heap_pos {new int[_size]},
  _h {h},
  _height {height},
  _width {width}
{
  init();
}

template <class base, int num_nb>
marcher<base, num_nb>::marcher(
  int height, int width, double h, double const * s_cache):
  _size {width*height},
  _heap {{this}, initial_heap_capacity(_size)},
  _U {new double[_size]},
  _s_cache {new double[_size]},
  _state {new state[_size]},
  _heap_pos {new int[_size]},
  _h {h},
  _height {height},
  _width {width}
{
  init();

  memcpy((void *) _s_cache, (void *) s_cache, sizeof(double)*height*width);
}

template <class base, int num_nb>
marcher<base, num_nb>::marcher(
  int height, int width, double h,
  std::function<double(double, double)> s, double x0, double y0):
  _size {width*height},
  _heap {{this}, initial_heap_capacity(_size)},
  _U {new double[_size]},
  _s_cache {new double[_size]},
  _state {new state[_size]},
  _heap_pos {new int[_size]},
  _h {h},
  _height {height},
  _width {width}
{
  init();

  double * ptr = const_cast<double *>(_s_cache);
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      ptr[linear_index(i, j)] = s(h*j - x0, h*i - y0);
    }
  }
}

template <class base, int num_nb>
marcher<base, num_nb>::~marcher()
{
  delete[] _U;
  delete[] _s_cache;
  delete[] _state;
  delete[] _heap_pos;
}

template <class base, int num_nb>
void marcher<base, num_nb>::init()
{
  for (int i = 0; i < _size; ++i) {
    _U[i] = inf<double>;
  }

  for (int i = 0; i < _size; ++i) {
    _state[i] = state::far;
  }
}

template <class base, int num_nb>
void marcher<base, num_nb>::run()
{
  while (!_heap.empty()) {
    int lin = _heap.front();
    _heap.pop_front();
    _state[lin] = state::valid;
    visit_neighbors(lin);
  }  
}

template <class base, int num_nb>
void
marcher<base, num_nb>::add_boundary_node(int i, int j, double value)
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
#endif
  int lin = linear_index(i, j);
  _U[lin] = value;
  _state[lin] = state::trial;
  _heap.insert(lin);
}

template <class base, int num_nb>
void
marcher<base, num_nb>::add_boundary_nodes(
  int const * i, int const * j, double const * U, int num)
{
  for (int k = 0; k < num; ++k) {
    add_boundary_node(i[k], j[k], U[k]);
  }
}

template <class base, int num_nb>
void
marcher<base, num_nb>::add_boundary_nodes(
  std::tuple<int, int, double> const * nodes, int num)
{
  for (int k = 0; k < num; ++k) {
    add_boundary_node(
      std::get<0>(nodes[k]),
      std::get<1>(nodes[k]),
      std::get<2>(nodes[k]));
  }
}

#define LINE(p0, u0, s, s0, h)                          \
  updates::line<base::F_>()(norm2<2>(p0), u0, s, s0, h)

template <class base, int num_nb>
void
marcher<base, num_nb>::add_boundary_node(
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

    int lin = linear_index(i_, j_);
    _U[lin] = LINE(P[a], u0, _s_cache[lin], s0, h);
    _state[lin] = state::trial;
    _heap.insert(lin);
  }
}

#undef LINE

template <class base, int num_nb>
void
marcher<base, num_nb>::set_fac_src(int i, int j, fac_src const * fc)
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
  assert(in_bounds(fc->i, fc->j));
#endif
  _lin2fac[linear_index(i, j)] = fc;
}

template <class base, int num_nb>
double
marcher<base, num_nb>::get_value(int i, int j) const
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
#endif
  return _U[linear_index(i, j)];
}

template <class base, int num_nb>
bool
marcher<base, num_nb>::in_bounds(int i, int j) const
{
  return (unsigned) i < (unsigned) _height && (unsigned) j < (unsigned) _width;
}

template <class base, int num_nb>
bool marcher<base, num_nb>::in_bounds(double i, double j) const {
  return 0 <= i <= _height - 1 && 0 <= j <= _width - 1;
}

template <class base, int num_nb>
double
marcher<base, num_nb>::get_speed(int i, int j) const
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j));
  assert(_s_cache != nullptr);
#endif
  return _s_cache[_width*i + j];
}

// TODO: we want to delete this---right now, it's a bit muddled, since
// it should probably be called `is_inbounds_and_valid'... but we want
// to get rid of `in_bounds', anyway, so once we do that, this
// function will just be deleted
template <class base, int num_nb>
bool
marcher<base, num_nb>::is_valid(int i, int j) const
{
  return in_bounds(i, j) && _state[linear_index(i, j)] == state::valid;
}

template <class base, int num_nb>
void
marcher<base, num_nb>::visit_neighbors(int lin_center)
{
  int const i = get_i(lin_center);
  int const j = get_j(lin_center);

  // These are temporary indices used below, analogous to i and j,
  // respectively.
  int a, b, lin;

  // Traverse the update neighborhood of n and set all far nodes to
  // trial and insert them into the heap.
  for (int k = 0; k < num_nb; ++k) {
    a = i + __di(k), b = j + __dj(k), lin = linear_index(a, b);
    if (in_bounds(a, b) && _state[lin] == state::far) {
      _state[lin] = state::trial;
      _heap.insert(lin);
    }
  }

  // Find the valid neighbors in the "full" neighborhood of n
  // (i.e. the unit max norm ball).
  for (int k = 0; k < 8; ++k) {
    valid_nb[k] = -1;
    a = i + __di(k), b = j + __dj(k), lin = linear_index(a, b);
    if (in_bounds(a, b) && _state[lin] == state::valid) {
      valid_nb[k] = lin;
    }
  }

  // Some explanation of the indices used below:
  // - l is a radial index circling (a, b)
  // - parent is the radial index of (i, j) expressed in the same
  //   index space as l
  int di_k, dj_k;
  int * nb = static_cast<base *>(this)->nb;
  auto const set_nb = [&] (int parent) {
    for (int l = 0; l < num_nb; ++l) {
      nb[l] = -1;
    }
    nb[parent] = lin_center;
    for (int l = 0; l < num_nb; ++l) {
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
        nb[l] = valid_nb[m];
      }
    }
  };

  // Update the node at (i, j). Before calling, `nb' needs to be
  // filled appropriately. Upon updating, this sets the value of n and
  // adjusts its position in the heap.
  auto const update = [&] (int lin_hat) {
    auto T = inf<double>;
    static_cast<base *>(this)->s_hat = _s_cache[lin_hat];
    update_impl(lin_hat, T);
    if (T < _U[lin_hat]) {
#if EIKONAL_DEBUG && !RELWITHDEBINFO
      assert(T >= 0);
#endif
      _U[lin_hat] = T;
      _heap.update(lin_hat);
    }
  };

  // Get the parent index of a radial index `k'.
  auto const get_parent = [] (int k) {
    if (num_nb == 4) {
      return (k + 2) % 4;
    } else {
      if (k < 4) return (k + 2) % 4;
      else return ((k - 2) % 4) + 4;
    }
  };

  // This is the main update loop. Each neighbor of n which isn't
  // valid is now trial. For each neighboring trial node, use
  // `set_nb' to grab its valid neighbors and use `update' to
  // actually update the node's value and update its position in the
  // heap.
  for (int k = 0; k < num_nb; ++k) {
    if (valid_nb[k] == -1) {
      di_k = __di(k), dj_k = __dj(k);
      a = i + di_k, b = j + dj_k;
      if (!in_bounds(a, b)) continue;
      int parent = get_parent(k);
      set_nb(parent);
      update(linear_index(a, b));
    }
  }
}

#undef __di
#undef __dj

#endif // __MARCHER_IMPL_HPP_HPP__
