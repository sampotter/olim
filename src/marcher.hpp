#pragma once

#include <src/config.hpp>

// TODO: remove these
#include <unordered_map>

#include "heap.hpp"
#include "slow.hpp"
#include "state.hpp"
#include "vec.hpp"

constexpr int max_num_nb(int n) {
  int lut[2] = {8, 26};
  return lut[n - 2];
}

template <int n>
struct fac_src
{
  fac_src(vec<double, n> coords, double s):
    coords {coords + decltype(coords)::one()}, s {s} {}

  vec<double, n> coords;
  double s;
};

template <class base, int n, int num_nb>
struct marcher
{
  using fac_src_t = fac_src<n>;

  // These are for use with our pybind11 bindings. They aren't used
  // internally.
  using float_type = double;

  using fvec = vec<double, n>;
  using ivec = vec<int, n>;
  using uvec = vec<unsigned, n>;

  static constexpr int ndim = n;

  marcher(ivec dims, double h, no_slow_t const &);
  marcher(ivec dims, double h);
  marcher(ivec dims, double h, double const * s);
  virtual ~marcher();

  void run();

  void add_boundary_node(int * inds, double U = 0.0);
  void add_boundary_node(ivec inds, double U = 0.0);
  void add_boundary_nodes(ivec const * inds, double const * U, int num);
  void add_boundary_node(fvec coords, double s, double U = 0.0);

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
    return ::to_linear_index(inds, _dims);
  }

  inline ivec to_vector_index(int lin) const {
    return ::to_vector_index(lin, _dims);
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
      return _m->_heap_pos[lin];
    }

    inline void set_heap_pos(int lin, int pos) {
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

  int _linear_offsets[max_num_nb(n)];

  // TODO: this is a quick hack just to get this working for the time
  // being.
  std::unordered_map<int, fac_src<n> const *> _lin2fac;
};

#include "marcher.impl.hpp"
