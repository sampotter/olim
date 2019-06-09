#pragma once

#include "heap.hpp"
#include "state.hpp"
#include "vec.hpp"

template <class T, int n, ordering ord>
void fill_boundary(vec<int, n> dims, T * ptr, T value)
{
  // First, iterate each of each dimension `d':
  for (int d = 0; d < n; ++d) {
    // Extract the dimensions not equal to `d' as a new dimension
    // vector called `subdims'.
    vec<int, n - 1> subdims;
    int j = 0;
    for (int i = 0; i < n; ++i) {
      if (i != d) {
        subdims[j] = dims[i];
        ++j;
      }
    }

    // Iterate over each point spanned by `subdims'.
    decltype(dims) inds;
    int lin;
    for (int i = 0; i < subdims.product(); ++i) {
      // Prepare an index vector of full dimension with indices
      // corresponding to `subdims'.
      vec<int, n - 1> subinds = to_vector_index<ord>(i, subdims);
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
      lin = to_linear_index<ord>(inds, dims);
      ptr[lin] = value;

      inds[d] = dims[d] - 1;
      lin = to_linear_index<ord>(inds, dims);
      ptr[lin] = value;
    }
  }
}

template <class base, int n>
struct base_marcher
{
  // TODO: duplicated...
  using ivec = vec<int, n>;
  using uvec = vec<unsigned, n>;

  base_marcher(ivec dims):
    _dims {dims + 2*ivec::one()},
    _size {_dims.product()},
    _heap {{this}, initial_heap_capacity(_size)},
    _U {new double[_size]},
    _state {new state[_size]},
    _heap_pos {new int[_size]}
  {
    for (int i = 0; i < _size; ++i) {
      _U[i] = inf<double>;
    }

    for (int i = 0; i < _size; ++i) {
      _state[i] = state::far;
    }
    fill_boundary<state, n, base::get_ord()>(_dims, _state, state::boundary);

#if OLIM_DEBUG && !RELWITHDEBINFO
    /**
     * If we're building in "slow debug mode", then we additionally
     * initialize the heap positions to -1 for use with debug asserts
     * elsewhere in `marcher'.
     */
    for (int i = 0; i < _size; ++i) {
      _heap_pos[i] = -1;
    }
#endif
  }

  ~base_marcher() {
    delete[] _U;
    delete[] _state;
    delete[] _heap_pos;
  }

  void solve();
  int step();
  int step_impl();
  bool peek(double * U, int * lin) const;
  void adjust(int const * inds, double U);
  void adjust(ivec inds, double U);

  double get_U(ivec inds) const;
  state get_state(ivec inds) const;

  inline double * get_U_ptr() const {
    return this->_U;
  }

  inline char * get_state_ptr() const {
    return reinterpret_cast<char *>(this->_state);
  }


OLIM_PROTECTED:

  bool in_bounds(ivec inds) const;

  struct proxy
  {
    using value_t = double;

    proxy(base_marcher * m): _m {m} {}

    inline value_t get_value(int lin) const {
      return _m->_U[lin];
    }

    inline int get_heap_pos(int lin) const {
#if OLIM_DEBUG && !RELWITHDEBINFO
      assert(_m->_heap_pos[lin] != -1);
#endif
      return _m->_heap_pos[lin];
    }

    inline void set_heap_pos(int lin, int pos) {
#if OLIM_DEBUG && !RELWITHDEBINFO
      assert(pos >= 0);
      assert(lin >= 0);
      assert(lin < _m->_size);
#endif
      _m->_heap_pos[lin] = pos;
    }

    base_marcher * _m {nullptr};
  };

  ivec _dims;
  int _size;

  heap<int, proxy> _heap;

  double * _U {nullptr};
  state * _state {nullptr};
  int * _heap_pos {nullptr};
};

#include "base_marcher.impl.hpp"
