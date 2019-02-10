#pragma once

#include <src/config.hpp>

// TODO: remove these
#include <unordered_map>

#include "heap.hpp"
#include "slow.hpp"
#include "state.hpp"
#include "vec.hpp"

template <int N>
struct fac_src
{
  fac_src(vec<double, N> coords, double s): coords {coords}, s {s} {}

  vec<double, N> coords;
  double s;
};

template <class base, int N, int num_nb>
struct marcher
{
  using fac_src_t = fac_src<N>;

  // These are for use with our pybind11 bindings. They aren't used
  // internally.
  using float_type = double;

  using fvec = vec<double, N>;
  using ivec = vec<int, N>;
  using uvec = vec<unsigned, N>;

  static constexpr int ndim = N;

  marcher(ivec dims, double h);
  marcher(ivec dims, double h, no_slow_t const &);
  marcher(ivec dims, double h, double const * s);
  marcher(ivec dims, double h, slow<N> s, fvec origin = fvec::zero());
  virtual ~marcher();

  void run();

  void add_boundary_node(ivec inds, double U = 0.0);
  void add_boundary_nodes(ivec const * inds, double const * U, int num);
  void add_boundary_node(fvec coords, double s, double U = 0.0);

  void set_fac_src(ivec inds, fac_src<N> const * src);

  double get_s(ivec inds) const;
  double get_U(ivec inds) const;

  inline double * get_U_ptr() const {
    printf("get_U_ptr(): %p\n", _U);
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
  bool is_valid(ivec inds) const;

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

  int _size;

  heap<int, proxy> _heap;
  
  double * _U {nullptr};
  double * _s {nullptr};
  state * _state {nullptr};
  int * _heap_pos {nullptr};
  double _h {1};
  ivec _dims;

  // TODO: this is a quick hack just to get this working for the time
  // being.
  std::unordered_map<int, fac_src<N> const *> _lin2fac;
};

#include "marcher.impl.hpp"
