#pragma once

#include <src/config.hpp>

// TODO: remove these
#include <unordered_map>

#include "fac.hpp"
#include "heap.hpp"
#include "slow.hpp"
#include "state.hpp"
#include "vec.hpp"

constexpr int max_num_nb(int n) {
  // TODO: would be nice to rewrite this as 3^n - 1, but need a
  // constexpr int power function
  int lut[2] = {8, 26};
  return lut[n - 2];
}

template <class base, int n, int num_nb, ordering ord = ordering::COLUMN_MAJOR>
struct marcher
{
  using fac_src_t = fac_src<n>;

  using fvec = vec<double, n>;
  using ivec = vec<int, n>;
  using uvec = vec<unsigned, n>;

  static constexpr int ndim = n;

  marcher(ivec dims, double h, no_slow_t const &);
  marcher(ivec dims, double h);
  marcher(ivec dims, double h, double const * s);
  virtual ~marcher();

  void solve();
  int step();
  int step_impl();

  inline double peek() const {
    return _U[_heap.front()];
  }

  void add_src(int * inds, double U = 0.0);
  void add_src(ivec inds, double U = 0.0);
  void add_srcs(ivec const * inds, double const * U, int num);
  void add_src(fvec coords, double s, double U = 0.0);

  void add_bd(int const * inds);
  void add_bd(ivec inds);

  void add_free(int const * inds);
  void add_free(ivec inds);

  void set_fac_src(ivec inds, fac_src<n> const * src);

  double get_U(ivec inds) const;
  double get_s(ivec inds) const;
  state get_state(ivec inds) const;

  inline double * get_U_ptr() const {
    return _U;
  }

  inline double * get_s_ptr() const {
    return _s;
  }

  inline char * get_state_ptr() const {
    return reinterpret_cast<char *>(_state);
  }

  inline void set_s_ptr(double * s) {
    _s = s;
  }

OLIM_PROTECTED:
  inline int to_linear_index(ivec inds) const {
    return ::to_linear_index<ord>(inds, _dims);
  }

  inline ivec to_vector_index(int lin) const {
    return ::to_vector_index<ord>(lin, _dims);
  }

  bool in_bounds(ivec inds) const;

  inline bool is_factored(int lin) const {
    return _lin2fac.find(lin) != _lin2fac.end();
  }

  double get_h() const { return _h; }
  
  void visit_neighbors(int lin);
  virtual void update_impl(int lin, int const * nb, int parent, double & U) = 0;

  struct proxy
  {
    using value_t = double;

    proxy(marcher * m): _m {m} {}

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

    marcher * _m {nullptr};
  };

  ivec _dims;
  int _size;

  heap<int, proxy> _heap;
  
  double * _U {nullptr};
  double * _s {nullptr};
  state * _state {nullptr};
  int * _heap_pos {nullptr};
  double _h {1};

  int _linear_offset[max_num_nb(n)];
  int _child_offset[num_nb][num_nb];

  // TODO: a quick hack just to get this working for now
  std::unordered_map<int, fac_src<n> const *> _lin2fac;
};

#include "marcher.impl.hpp"
