#pragma once

#include <assert.h>
#include <math.h>

#include "common.hpp"
#include "offsets.hpp"
#include "range.hpp"
#include "updates.line.hpp"

// TODO: really need an external memory constructor that will let us
// use external memory somewhere for _s and _U so we don't have
// to double up on these...

template <class base, int n, int num_nb>
marcher<base, n, num_nb>::marcher(ivec dims, double h, no_slow_t const &):
  _size {dims.product()},
  _heap {{this}, initial_heap_capacity(_size)},
  _U {new double[_size]},
  _s {new double[_size]},
  _state {new state[_size]},
  _heap_pos {new int[_size]},
  _h {h},
  _dims {dims}
{
  for (int i = 0; i < _size; ++i) {
    _U[i] = inf<double>;
  }

  for (int i = 0; i < _size; ++i) {
    _state[i] = state::far;
  }
}

template <class base, int n, int num_nb>
marcher<base, n, num_nb>::marcher(ivec dims, double h):
  marcher {dims, h, no_slow_t {}}
{
  for (int i = 0; i < _size; ++i) {
    _s[i] = 1;
  }
}

template <class base, int n, int num_nb>
marcher<base, n, num_nb>::marcher(ivec dims, double h, double const * s):
  marcher {dims, h, no_slow_t {}}
{
  memcpy((void *) _s, (void *) s, sizeof(double)*_size);
}

template <class base, int n, int num_nb>
marcher<base, n, num_nb>::marcher(ivec dims, double h, slow<n> s, fvec origin):
  marcher {dims, h, no_slow_t {}}
{
  for (auto inds: range<n> {dims}) {
    _s[to_linear_index(inds)] = s(h*inds - origin);
  }
}

template <class base, int n, int num_nb>
marcher<base, n, num_nb>::~marcher()
{
  delete[] _U;
  delete[] _s;
  delete[] _state;
  delete[] _heap_pos;
}

template <class base, int n, int num_nb>
void marcher<base, n, num_nb>::run()
{
  while (!_heap.empty()) {
    int lin = _heap.front();
    _heap.pop_front();
    _state[lin] = state::valid;
    visit_neighbors(lin);
  }  
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_boundary_node(int * inds, double U)
{
  add_boundary_node(ivec {inds}, U);
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_boundary_node(ivec inds, double U)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  int lin = to_linear_index(inds);
  _U[lin] = U;
  _state[lin] = state::trial;
  _heap.insert(lin);
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_boundary_nodes(ivec const * inds, double const * U, int num)
{
  for (int i = 0; i < num; ++i) {
    add_boundary_node(inds[i], U[i]);
  }
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_boundary_node(fvec coords, double s, double U)
{
  double h = get_h(), u0 = U, s0 = s;
  fvec inds = coords/h;
  assert(in_bounds(inds));

  // TODO: this isn't as general as it could be. We also want to
  // handle cases where i or j are grid-aligned, in which case we need
  // to process 6 nodes. If i and j are both grid aligned, then we
  // should just call the integer add_boundary_node.
  int corners[n][2];
  for (int i = 0; i < n; ++i) {
    corners[i][0] = (int) floor(inds[i]);
    corners[i][1] = (int) floor(inds[i]) + 1;
  }

  for (auto inds_: range<n>::cube(2)) {
    ivec inds__;
    for (int i = 0; i < n; ++i) {
      inds__[i] = corners[i][inds_[i]];
    }
    assert(in_bounds(inds__));

    fvec p = fvec {inds} - inds__;
    int lin = to_linear_index(inds__);
    _U[lin] = updates::line<base::F_>()(p.norm2(), u0, _s[lin], s0, h);
    _state[lin] = state::trial;
    _heap.insert(lin);
  }
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::set_fac_src(ivec inds, fac_src<n> const * fc)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  _lin2fac[to_linear_index(inds)] = fc;
}

template <class base, int n, int num_nb>
double
marcher<base, n, num_nb>::get_U(ivec inds) const
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  return _U[to_linear_index(inds)];
}

template <class base, int n, int num_nb>
bool
marcher<base, n, num_nb>::in_bounds(ivec inds) const
{
  return uvec {inds} < uvec {_dims};
}

template <class base, int n, int num_nb>
double
marcher<base, n, num_nb>::get_s(ivec inds) const
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
  assert(_s != nullptr);
#endif
  return _s[to_linear_index(inds)];
}

// TODO: we want to delete this---right now, it's a bit muddled, since
// it should probably be called `is_inbounds_and_valid'... but we want
// to get rid of `in_bounds', anyway, so once we do that, this
// function will just be deleted
template <class base, int n, int num_nb>
bool
marcher<base, n, num_nb>::is_valid(ivec inds) const
{
  return in_bounds(inds) && _state[to_linear_index(inds)] == state::valid;
}

constexpr int max_num_nb(int n) {
  int lut[2] = {8, 26};
  return lut[n - 2];
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::visit_neighbors(int lin_center)
{
  ivec inds = to_vector_index(lin_center);

  int valid_nb[max_num_nb(n)];
  int child_nb[num_nb];

  // Traverse the update neighborhood of n and set all far nodes to
  // trial and insert them into the heap.
  for (int i = 0; i < num_nb; ++i) {
    ivec inds_ = inds + get_offset<n>(i);
    if (in_bounds(inds_)) {
      int lin = to_linear_index(inds_);
      if (_state[lin] == state::far) {
        _state[lin] = state::trial;
        _heap.insert(lin);
      }
    }
  }

  // Find the valid neighbors in the "full" neighborhood of n
  // (i.e. the unit max norm ball).
  for (int i = 0; i < max_num_nb(n); ++i) {
    valid_nb[i] = -1;
    ivec inds_ = inds + get_offset<n>(i);
    if (in_bounds(inds_))  {
      int lin = to_linear_index(inds_);
      if (_state[lin] == state::valid) {
        valid_nb[i] = lin;
      }
    }
  }

  // Some explanation of the indices used below:
  // - l is a radial index circling (a, b)
  // - parent is the radial index of (i, j) expressed in the same
  //   index space as l
  auto const set_child_nb = [&] (int parent, ivec offset) {
    for (int l = 0; l < num_nb; ++l) {
      child_nb[l] = -1;
    }
    child_nb[parent] = lin_center;
    for (int l = 0; l < num_nb; ++l) {
      if (l == parent) {
        continue;
      }
      ivec offset_ = offset + get_offset<n>(l);
      if (offset_.normi() > 1) {
        continue;
      }
      if (in_bounds(inds + offset_)) {
        child_nb[l] = valid_nb[get_linear_offset<n>(offset_)];
      }
    }
  };

  // Update the node at (i, j). Before calling, `nb' needs to be
  // filled appropriately. Upon updating, this sets the value of n and
  // adjusts its position in the heap.
  auto const update = [&] (int lin_hat, int parent) {
    auto U = inf<double>;
    update_impl(lin_hat, child_nb, parent, U);
    if (U < _U[lin_hat]) {
#if OLIM_DEBUG && !RELWITHDEBINFO
      assert(U >= 0);
#endif
      _U[lin_hat] = U;
      _heap.update(lin_hat);
    }
  };

  // Get the parent index of a radial index `i'.
  //
  // TODO: `get_parent' is a contender for the worst code in the whole
  // project. Would it be better to just use a LUT here? Should try
  // this at some point...
  auto const get_parent = [] (int i) {
    if (n == 2) {
      if (num_nb <= 4) {
        return (i + 2) % 4;
      } else {
        if (i < 4) return (i + 2) % 4;
        else return ((i - 2) % 4) + 4;
      }
    } else if (n == 3) {
      if (num_nb <= 6) {
        return (i + 3) % 6;
      } else if (num_nb <= 18) {
        return i < 6 ? (i + 3) % 6 : 22 - 2*(i/2) + (i % 2);
      } else {
        if (i < 6) return (i + 3) % 6;
        else if (i < 18) return 22 - 2*(i/2) + (i % 2);
        else return 42 - 2*(i/2) + (i % 2);
      }
    }
  };

  // This is the main update loop. Each neighbor of n which isn't
  // valid is now trial. For each neighboring trial node, use
  // `set_child_nb' to grab its valid neighbors and use `update' to
  // actually update the node's value and update its position in the
  // heap.
  for (int i = 0; i < num_nb; ++i) {
    if (valid_nb[i] == -1) {
      ivec offset = get_offset<n>(i);
      ivec inds_ = inds + offset;
      if (in_bounds(inds_)) {
        int parent = get_parent(i);
        set_child_nb(parent, offset);
        update(to_linear_index(inds_), parent);
      }
    }
  }
}
