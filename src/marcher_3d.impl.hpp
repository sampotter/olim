#pragma once

#include <src/config.hpp>

#include <assert.h>

#include "offsets.hpp"
#include "updates.line.hpp"

template <class base, int num_nb>
marcher_3d<base, num_nb>::marcher_3d() {}

template <class base, int num_nb>
marcher_3d<base, num_nb>::marcher_3d(vec3<int> dims, double h, no_slow_t const &):
  _size {dims.product()},
  _heap {{this}, initial_heap_capacity(_size)},
  _U {new double[_size]},
  _s_cache {new double[_size]},
  _state {new state[_size]},
  _heap_pos {new int[_size]},
  _h {h},
  _dims {dims}
{
  init();
}

template <class base, int num_nb>
marcher_3d<base, num_nb>::marcher_3d(vec3<int> dims, double h,
                                     std::function<double(double, double, double)> s,
                                     double x0, double y0, double z0):
  _size {dims.product()},
  _heap {{this}, initial_heap_capacity(_size)},
  _U {new double[_size]},
  _s_cache {new double[_size]},
  _state {new state[_size]},
  _heap_pos {new int[_size]},
  _h {h},
  _dims {dims}
{
  init();

  // Grab a writable pointer to cache the values of `s'.
  double x, z, * ptr = const_cast<double *>(_s_cache);
  for (int k = 0; k < _dims[2]; ++k) {
    z = h*k - z0;
    for (int j = 0; j < _dims[1]; ++j) {
      x = h*j - x0;
      for (int i = 0; i < _dims[0]; ++i) {
        ptr[linear_index(i, j, k)] = s(x, h*i - y0, z);
      }
    }
  }
}

template <class base, int num_nb>
marcher_3d<base, num_nb>::marcher_3d(vec3<int> dims, double h, double const * s_cache):
  _size {dims.product()},
  _heap {{this}, initial_heap_capacity(_size)},
  _U {new double[_size]},
  _s_cache {new double[_size]},
  _state {new state[_size]},
  _heap_pos {new int[_size]},
  _h {h},
  _dims {dims}
{
  init();

  memcpy((void *) _s_cache, (void *) s_cache, sizeof(double)*_size);
}

template <class base, int num_nb>
marcher_3d<base, num_nb>::~marcher_3d()
{
  delete[] _U;
  delete[] _s_cache;
  delete[] _state;
  delete[] _heap_pos;
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::init()
{
  for (int i = 0; i < _size; ++i) {
    _U[i] = inf<double>;
  }

  for (int i = 0; i < _size; ++i) {
    _state[i] = state::far;
  }
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::run()
{
  while (!_heap.empty()) {
    int lin = _heap.front();
    _heap.pop_front();
    _state[lin] = state::valid;
    visit_neighbors(lin);
  }  
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::add_boundary_node(
  int i, int j, int k, double value)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
#endif
  int lin = linear_index(i, j, k);
  _U[lin] = value;
  _state[lin] = state::trial;
  _heap.insert(lin);
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::add_boundary_nodes(
  int const * i, int const * j, int const * k, double const * U, int num)
{
  for (int l = 0; l < num; ++l) {
    add_boundary_node(i[l], j[l], k[l], U[l]);
  }
}

template <class base, int num_nb>
void
marcher_3d<base, num_nb>::add_boundary_node(
  double x, double y, double z, double s, double value)
{
  double h = get_h(), i = y/h, j = x/h, k = z/h, u0 = value, s0 = s;
  assert(in_bounds(i, j, k));

  // TODO: make this more general: see comment in marcher.impl.hpp for
  // related function
  int is[2] = {(int) floor(i), (int) floor(i) + 1};
  int js[2] = {(int) floor(j), (int) floor(j) + 1};
  int ks[2] = {(int) floor(k), (int) floor(k) + 1};

  vec3<double> P[8] = {
    {i - is[0], j - js[0], k - ks[0]},
    {i - is[1], j - js[0], k - ks[0]},
    {i - is[0], j - js[1], k - ks[0]},
    {i - is[1], j - js[1], k - ks[0]},
    {i - is[0], j - js[0], k - ks[1]},
    {i - is[1], j - js[0], k - ks[1]},
    {i - is[0], j - js[1], k - ks[1]},
    {i - is[1], j - js[1], k - ks[1]},
  };

  for (int a = 0; a < 8; ++a) {
    int b0 = a & 1, b1 = (a & 2) >> 1, b2 = (a & 4) >> 2;
    int i_ = is[b0], j_ = js[b1], k_ = ks[b2];

    assert(in_bounds(i_, j_, k_));

    int lin = linear_index(i_, j_, k_);
    _U[lin] = updates::line<base::F_>()(P[a].norm2(), u0, _s_cache[lin], s0, h);
    _state[lin] = state::trial;
    _heap.insert(lin);
  }
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::set_fac_src(int i, int j, int k, fac_src_3d const * fc)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
  assert(in_bounds(fc->i, fc->j, fc->k));
#endif
  _lin2fac[linear_index(i, j, k)] = fc;
}

template <class base, int num_nb>
double marcher_3d<base, num_nb>::get_value(int i, int j, int k) const {
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
#endif
  return _U[linear_index(i, j, k)];
}

template <class base, int num_nb>
bool marcher_3d<base, num_nb>::in_bounds(int i, int j, int k) const {
  return (unsigned) i < (unsigned) _dims[0] &&
    (unsigned) j < (unsigned) _dims[1] && (unsigned) k < (unsigned) _dims[2];
}

template <class base, int num_nb>
double marcher_3d<base, num_nb>::get_s(int i, int j, int k) const {
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(i, j, k));
  assert(_s_cache != nullptr);
#endif
  return _s_cache[linear_index(i, j, k)];
}

// TODO: we want to delete this---right now, it's a bit muddled, since
// it should probably be called `is_inbounds_and_valid'... but we want
// to get rid of `in_bounds', anyway, so once we do that, this
// function will just be deleted
template <class base, int num_nb>
bool marcher_3d<base, num_nb>::is_valid(int i, int j, int k) const {
  return in_bounds(i, j, k) && _state[linear_index(i, k, k)] == state::valid;
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::visit_neighbors(int lin_center) {
  int const i = get_i(lin_center);
  int const j = get_j(lin_center);
  int const k = get_k(lin_center);

  // See comments in marcher.impl.hpp; the visit_neighbors_impl there
  // is done analogously to this one.

  int a, b, c, lin;

  // Stage neighbors.
  for (int l = 0; l < num_nb; ++l) {
    a = i + di<3>[l], b = j + dj<3>[l], c = k + dk<3>[l],
      lin = linear_index(a, b, c);
    if (in_bounds(a, b, c) && _state[lin] == state::far) {
      _state[lin] = state::trial;
      _heap.insert(lin);
    }
  }

  // Get valid neighbors.
  int valid_nb[26], child_nb[num_nb];
  for (int l = 0; l < 26; ++l) {
    valid_nb[l] = -1;
    a = i + di<3>[l], b = j + dj<3>[l], c = k + dk<3>[l];
    lin = linear_index(a, b, c);
    if (in_bounds(a, b, c) && _state[lin] == state::valid) {
      valid_nb[l] = lin;
    }
  }

  int di_l, dj_l, dk_l;
  auto const set_child_nb = [&] (int parent) {
    for (int m = 0; m < num_nb; ++m) {
      child_nb[m] = -1;
    }
    child_nb[parent] = lin_center;
    for (int m = 0; m < num_nb; ++m) {
      if (m == parent) {
        continue;
      }
      int di_lm, dj_lm, dk_lm;
      if (std::abs(di_lm = di_l + di<3>[m]) > 1 ||
          std::abs(dj_lm = dj_l + dj<3>[m]) > 1 ||
          std::abs(dk_lm = dk_l + dk<3>[m]) > 1) {
        continue;
      }
      if (in_bounds(i + di_lm, j + dj_lm, k + dk_lm)) {
        int p = d2l(di_lm, dj_lm, dk_lm);
        child_nb[m] = valid_nb[p];
      }
    }
  };

  auto & s_hat = static_cast<base *>(this)->s_hat;
  auto const update = [&] (int lin_hat, int parent) {
    auto U = inf<double>;
    s_hat = _s_cache[lin_hat];
    update_impl(lin_hat, child_nb, parent, U);
    if (U < _U[lin_hat]) {
#if OLIM_DEBUG && !RELWITHDEBINFO
      assert(U >= 0);
#endif
      _U[lin_hat] = U;
      _heap.update(lin_hat);
    }
  };

  auto const get_parent = [] (int l) {
    // TODO: check base::num_nb to reduce amount of branching here
    if (l < 6) return (l + 3) % 6;
    else if (l < 18) return 22 - 2*(l/2) + (l % 2);
    else return 42 - 2*(l/2) + (l % 2);
  };

  for (int l = 0; l < num_nb; ++l) {
    if (valid_nb[l] == -1) {
      di_l = di<3>[l], dj_l = dj<3>[l], dk_l = dk<3>[l];
      a = i + di_l, b = j + dj_l, c = k + dk_l;
      if (!in_bounds(a, b, c)) continue;
      int parent = get_parent(l);
      set_child_nb(parent);
      update(linear_index(a, b, c), parent);
    }
  }
}
