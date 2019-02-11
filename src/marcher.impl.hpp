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
  _dims {dims + 2*ivec::one()},
  _size {_dims.product()},
  _heap {{this}, initial_heap_capacity(_size)},
  _U {new double[_size]},
  _s {new double[_size]},
  _state {new state[_size]},
  _heap_pos {new int[_size]},
  _h {h}
{
  for (int i = 0; i < _size; ++i) {
    _U[i] = inf<double>;
  }

  for (int i = 0; i < _size; ++i) {
    _state[i] = state::far;
  }

  /**
   * Set the state of all nodes on the edge of the domain to
   * `boundary' and set boundary slowness values to infinity.
   *
   * First, iterate over each dimension `d':
   */
  for (int d = 0; d < n; ++d) {
    // Extract the dimensions not equal to `d' as a new dimension
    // vector called `subdims'.
    vec<int, n - 1> subdims;
    int j = 0;
    for (int i = 0; i < n; ++i) {
      if (i != d) {
        subdims[j] = _dims[i];
        ++j;
      }
    }

    // Iterate over each point spanned by `subdims'.
    ivec inds;
    int lin;
    for (int i = 0; i < subdims.product(); ++i) {
      // Prepare an index vector of full dimension with indices
      // corresponding to `subdims'.
      vec<int, n - 1> subinds = ::to_vector_index(i, subdims);
      int k = 0;
      for (int j = 0; j < n; ++j) {
        if (j != d) {
          inds[j] = subinds[k++];
        }
      }

      // Set dimension `d' to its two extremal values and set the
      // corresponding points' states to `boundary' and slowness
      // values to infinity.

      inds[d] = 0;
      lin = to_linear_index(inds);
      _state[lin] = state::boundary;
      _s[lin] = inf<double>;

      inds[d] = _dims[d] - 1;
      lin = to_linear_index(inds);
      _state[lin] = state::boundary;
      _s[lin] = inf<double>;
    }
  }

  /**
   * Precompute linear offsets. These are used in `visit_neighbors'.
   */
  for (int i = 0; i < max_num_nb(n); ++i) {
    _linear_offsets[i] = to_linear_index(get_offset<n>(i));
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
  // TODO: this isn't the most efficient way to do this. It would be
  // better to iterate over pointers to the heads of each maximum
  // contiguous block of memory and call memcpy for each one. This is
  // more complicated but would reduce the amount of index math
  // substantially.
  for (auto inds: range<n> {dims}) {
    _s[to_linear_index(inds + ivec::one())] = s[::to_linear_index(inds, dims)];
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
marcher<base, n, num_nb>::add_src(int * inds, double U)
{
  add_src(ivec {inds}, U);
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_src(ivec inds, double U)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  inds += ivec::one();
  int lin = to_linear_index(inds);
  _U[lin] = U;
  _state[lin] = state::trial;
  _heap.insert(lin);
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_srcs(ivec const * inds, double const * U, int num)
{
  for (int i = 0; i < num; ++i) {
    add_src(inds[i], U[i]);
  }
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_src(fvec coords, double s, double U)
{
  double h = get_h();
  fvec inds = coords/h;
  assert(in_bounds(inds));
  inds += fvec::one();

  // TODO: this isn't as general as it could be. We also want to
  // handle cases where i or j are grid-aligned, in which case we need
  // to process 6 nodes. If i and j are both grid aligned, then we
  // should just call the integer add_src.
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

    fvec p = inds - fvec {inds__};

    int lin = to_linear_index(inds__);
    _U[lin] = updates::line<base::F_>()(p.norm2(), U, _s[lin], s, h);
    _state[lin] = state::trial;
    _heap.insert(lin);
  }
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_bd(int * inds)
{
  add_bd(ivec {inds});
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::add_bd(ivec inds)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  inds += ivec::one();
  int lin = to_linear_index(inds);
  _U[lin] = inf<double>;
  _s[lin] = inf<double>;
  _state[lin] = state::boundary;
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::set_fac_src(ivec inds, fac_src<n> const * fc)
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(in_bounds(inds));
#endif
  inds += ivec::one();
  _lin2fac[to_linear_index(inds)] = fc;
}

template <class base, int n, int num_nb>
bool
marcher<base, n, num_nb>::in_bounds(ivec inds) const
{
  return uvec {inds} < uvec {_dims - 2*ivec::one()};
}

template <class base, int n, int num_nb>
double
marcher<base, n, num_nb>::get_U(ivec inds) const
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(_U != nullptr);
  assert(in_bounds(inds));
#endif
  return _U[to_linear_index(inds + ivec::one())];
}

template <class base, int n, int num_nb>
double
marcher<base, n, num_nb>::get_s(ivec inds) const
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(_s != nullptr);
  assert(in_bounds(inds));
#endif
  return _s[to_linear_index(inds + ivec::one())];
}

template <class base, int n, int num_nb>
state
marcher<base, n, num_nb>::get_state(ivec inds) const
{
#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(_state != nullptr);
  assert(in_bounds(inds));
#endif
  return _state[to_linear_index(inds + ivec::one())];
}

template <class base, int n, int num_nb>
void
marcher<base, n, num_nb>::visit_neighbors(int lin_center)
{
  int valid_nb[max_num_nb(n)];
  int child_nb[num_nb];

  // Traverse the update neighborhood of n and set all far nodes to
  // trial and insert them into the heap.
  for (int i = 0; i < num_nb; ++i) {
    int lin = lin_center + _linear_offsets[i];
    if (_state[lin] == state::far) {
      _state[lin] = state::trial;
      _heap.insert(lin);
    }
  }

  // Find the valid neighbors in the "full" neighborhood of n
  // (i.e. the unit max norm ball).
  for (int i = 0; i < max_num_nb(n); ++i) {
    valid_nb[i] = -1;
    int lin = lin_center + _linear_offsets[i];
    if (_state[lin] == state::valid) {
      valid_nb[i] = lin;
    }
  }

  // TODO: comment this
  auto const set_child_nb = [&] (int parent, ivec offset) {
    for (int i = 0; i < num_nb; ++i) {
      child_nb[i] = -1;
    }
    child_nb[parent] = lin_center;
    for (int i = 0; i < num_nb; ++i) {
      if (i == parent) {
        continue;
      }
      ivec offset_ = offset + get_offset<n>(i);
      if (offset_.normi() > 1) {
        continue;
      }
      child_nb[i] = valid_nb[get_linear_offset<n>(offset_)];
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
  // `valid' is now `trial' or `boundary'. We ignore `boundary' nodes;
  // for each neighboring trial node, use `set_child_nb' to grab its
  // valid neighbors and use `update' to actually update the node's
  // value and update its position in the heap.
  for (int i = 0; i < num_nb; ++i) {
    if (valid_nb[i] == -1) {
      int lin = lin_center + _linear_offsets[i];
      if (_state[lin] == state::trial) {
        int parent = get_parent(i);
        set_child_nb(parent, get_offset<n>(i));
        update(lin, parent);
      }
    }
  }
}
