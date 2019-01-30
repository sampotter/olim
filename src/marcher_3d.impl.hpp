#pragma once

#include <src/config.hpp>

#include <assert.h>

#include "offsets.hpp"
#include "range.hpp"
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
                                     std::function<double(vec3<double>)> s,
                                     vec3<double> origin):
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

  double * ptr = const_cast<double *>(_s_cache);
  for (auto inds: range<3> {_dims}) {
    ptr[to_linear_index(inds)] = s(h*inds - origin);
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
void marcher_3d<base, num_nb>::add_boundary_node(vec3<int> inds, double value)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  int lin = to_linear_index(inds);
  _U[lin] = value;
  _state[lin] = state::trial;
  _heap.insert(lin);
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::add_boundary_nodes(vec3<int> const * inds, double const * U, int num)
{
  for (int i = 0; i < num; ++i) {
    add_boundary_node(inds[i], U[i]);
  }
}

template <class base, int num_nb>
void
marcher_3d<base, num_nb>::add_boundary_node(vec3<double> coords, double s, double value)
{
  double h = get_h(), u0 = value, s0 = s;
  vec3<double> inds = coords/h;
  assert(in_bounds(inds));

  double i = inds[0], j = inds[1], k = inds[2];

  // TODO: make this more general: see comment in marcher.impl.hpp for
  // related function
  int is[2] = {(int) floor(i), (int) floor(i) + 1};
  int js[2] = {(int) floor(j), (int) floor(j) + 1};
  int ks[2] = {(int) floor(k), (int) floor(k) + 1};

  for (auto inds_: range<3> {{2, 2, 2}}) {
    vec3<int> inds__ = {is[inds_[0]], js[inds_[1]], ks[inds_[2]]};
    assert(in_bounds(inds__));

    vec3<double> p = vec3<double> {inds} - inds__;

    int lin = to_linear_index(inds__);
    _U[lin] = updates::line<base::F_>()(p.norm2(), u0, _s_cache[lin], s0, h);
    _state[lin] = state::trial;
    _heap.insert(lin);
  }
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::set_fac_src(vec3<int> inds, fac_src_3d const * fc)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  _lin2fac[to_linear_index(inds)] = fc;
}

template <class base, int num_nb>
double marcher_3d<base, num_nb>::get_value(vec3<int> inds) const {
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  return _U[to_linear_index(inds)];
}

template <class base, int num_nb>
bool marcher_3d<base, num_nb>::in_bounds(vec3<int> inds) const {
  return vec3<unsigned> {inds} < vec3<unsigned> {_dims};
}

template <class base, int num_nb>
double marcher_3d<base, num_nb>::get_s(vec3<int> inds) const {
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
  assert(_s_cache != nullptr);
#endif
  return _s_cache[to_linear_index(inds)];
}

// TODO: we want to delete this---right now, it's a bit muddled, since
// it should probably be called `is_inbounds_and_valid'... but we want
// to get rid of `in_bounds', anyway, so once we do that, this
// function will just be deleted
template <class base, int num_nb>
bool marcher_3d<base, num_nb>::is_valid(vec3<int> inds) const {
  return in_bounds(inds) && _state[to_linear_index(inds)] == state::valid;
}

template <class base, int num_nb>
void marcher_3d<base, num_nb>::visit_neighbors(int lin_center) {
  vec3<int> inds = to_vector_index(lin_center);

  // Stage neighbors.
  for (int i = 0; i < num_nb; ++i) {
    vec3<int> inds_ = inds + get_offset<3>(i);
    int lin = to_linear_index(inds_);
    if (in_bounds(inds_) && _state[lin] == state::far) {
      _state[lin] = state::trial;
      _heap.insert(lin);
    }
  }

  // Get valid neighbors.
  int valid_nb[26], child_nb[num_nb];
  for (int i = 0; i < 26; ++i) {
    valid_nb[i] = -1;
    vec3<int> inds_ = inds + get_offset<3>(i);
    int lin = to_linear_index(inds_);
    if (in_bounds(inds_) && _state[lin] == state::valid) {
      valid_nb[i] = lin;
    }
  }

  auto const set_child_nb = [&] (int parent, vec3<int> offset) {
    for (int m = 0; m < num_nb; ++m) {
      child_nb[m] = -1;
    }
    child_nb[parent] = lin_center;
    for (int m = 0; m < num_nb; ++m) {
      if (m == parent) {
        continue;
      }
      vec3<int> offset_ = offset + get_offset<3>(m);
      if (offset_.normi() > 1) {
        continue;
      }
      if (in_bounds(inds + offset_)) {
        child_nb[m] = valid_nb[get_linear_offset<3>(offset_)];
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

  for (int i = 0; i < num_nb; ++i) {
    if (valid_nb[i] == -1) {
      vec3<int> offset = get_offset<3>(i);
      vec3<int> inds_ = inds + offset;
      if (!in_bounds(inds_)) continue;
      int parent = get_parent(i);
      set_child_nb(parent, offset);
      update(to_linear_index(inds_), parent);
    }
  }
}
