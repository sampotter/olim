#pragma once

#include "heap.hpp"
#include "offsets.hpp"
#include "state.hpp"
#include "vec.hpp"

namespace detail {

constexpr int max_num_nb(int n) {
  // TODO: would be nice to rewrite this as 3^n - 1, but need a
  // constexpr int power function
  int lut[2] = {8, 26};
  return lut[n - 2];
}

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

}

template <class derived, int n>
struct base_marcher
{
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
    detail::fill_boundary<state, n, derived::get_ord()>(
      _dims, _state, state::boundary);

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

    /**
     * Precompute linear offsets. These are used in `visit_neighbors'.
     */
    for (int i = 0; i < detail::max_num_nb(n); ++i) {
      _linear_offset[i] = to_linear_index(detail::get_offset<n>(i));
    }
  }

  ~base_marcher() {
    delete[] _U;
    delete[] _state;
    delete[] _heap_pos;
  }

  void solve();
  int step();
  int step_impl();
  void add_src(int const * inds, double U = 0.0);
  void add_src(ivec inds, double U = 0.0);
  void add_srcs(ivec const * inds, double const * U, int num);
  void add_bd(int const * inds);
  void add_bd(ivec inds);
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


  void stage_neighbors(int lin);

  inline int to_linear_index(ivec inds) const {
    return ::to_linear_index<derived::get_ord()>(inds, this->_dims);
  }

  inline ivec to_vector_index(int lin) const {
    return ::to_vector_index<derived::get_ord()>(lin, this->_dims);
  }

  inline int to_external_linear_index(int lin) const {
    auto const inds = to_vector_index(lin) - ivec::one();
    return ::to_linear_index<derived::get_ord()>(inds, this->_dims - 2*ivec::one());
  }

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

  int _linear_offset[detail::max_num_nb(n)];
};

#include "base_marcher.impl.hpp"
