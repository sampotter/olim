#pragma once

#include <src/config.hpp>

// TODO: remove these
#include <functional>
#include <unordered_map>

#include "heap.hpp"
#include "slow.hpp"
#include "typedefs.h"
#include "vec.hpp"

struct fac_src
{
  fac_src(vec2<double> coords, double s): coords {coords}, s {s} {}

  vec2<double> coords;
  double s;
};

template <class base, int num_nb>
struct marcher
{
  using fac_src_t = fac_src;

  // These are for use with our pybind11 bindings. They aren't used
  // internally.
  using float_type = double;

  static constexpr int ndim = 2;
  // static constexpr int num_neighbors = base::num_neighbors;

  marcher(vec2<int> dims, double h, no_slow_t const &);
  marcher(vec2<int> dims, double h, double const * s_cache);
  marcher(vec2<int> dims, double h = 1,
          std::function<double(double, double)> s = static_cast<slow2>(s0),
          vec2<double> origin = vec2<double>::zero());
  virtual ~marcher();

  void init();

  void run();

  void add_boundary_node(vec2<int> inds, double value = 0.0);
  void add_boundary_nodes(vec2<int> const * inds, double const * U, int num);
  void add_boundary_nodes(std::tuple<int, int, double> const * nodes, int num);
  void add_boundary_node(vec2<double> coords, double s, double value = 0.0);

  void set_fac_src(vec2<int> inds, fac_src const * src);

  double get_s(vec2<int> inds) const;
  double get_value(vec2<int> inds) const;

  inline double const * get_s_cache_data() const {
    return _s_cache;
  }

OLIM_PROTECTED:
  inline int to_linear_index(vec2<int> inds) const {
    return ::to_linear_index(inds, _dims);
  }

  inline vec2<int> to_vector_index(int lin) const {
    return ::to_vector_index(lin, _dims);
  }

  bool in_bounds(vec2<int> inds) const;
  bool is_valid(vec2<int> inds) const;

  inline bool is_factored(int lin) const {
    return _lin2fac.find(lin) != _lin2fac.end();
  }

  double get_h() const { return _h; }
  
  void visit_neighbors(int lin);
  virtual void update_impl(int lin, double & U) = 0;

  int valid_nb[8];

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
  double const * _s_cache {nullptr};
  state * _state {nullptr};
  int * _heap_pos {nullptr};
  double _h {1};
  vec2<int> _dims;

  // TODO: this is a quick hack just to get this working for the time
  // being.
  std::unordered_map<int, fac_src const *> _lin2fac;
};

#include "marcher.impl.hpp"
