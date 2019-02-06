#pragma once

// TODO: try to remove this
#include <functional>
#include <unordered_map>

#include "heap.hpp"
#include "slow.hpp"
#include "typedefs.h"
#include "vec.hpp"

struct fac_src_3d
{
  fac_src_3d(vec3<double> coords, double s): coords {coords}, s {s} {}

  vec3<double> coords;
  double s;
};

template <class base, int num_nb>
struct marcher_3d
{
  using fac_src_t = fac_src_3d;

  // These are for use with our pybind11 bindings. They aren't used
  // internally.
  using float_type = double;

  static constexpr int ndim = 3;

  marcher_3d();
  marcher_3d(vec3<int> dims, double h, no_slow_t const &);
  marcher_3d(vec3<int> dims, double h = 1,
             std::function<double(vec3<double>)> s = static_cast<slow<3>>(s0<3>),
             vec3<double> origin = vec3<double>::zero());
  marcher_3d(vec3<int> dims, double h, double const * s);
  virtual ~marcher_3d();

  void init();

  void run();

  void add_boundary_node(vec3<int> inds, double value = 0.0);
  void add_boundary_nodes(vec3<int> const * inds, double const * U, int num);
  void add_boundary_node(vec3<double> coords, double s, double value = 0.0);

  void set_fac_src(vec3<int> inds, fac_src_3d const * src);

  double get_s(vec3<int> inds) const;
  double get_value(vec3<int> inds) const;

  inline double const * get_s_ptr() const {
    return _s;
  }

OLIM_PROTECTED:
  inline int to_linear_index(vec3<int> inds) const {
    return ::to_linear_index(inds, _dims);
  }

  inline vec3<int> to_vector_index(int lin) const {
    return ::to_vector_index(lin, _dims);
  }

  bool in_bounds(vec3<int> inds) const;
  bool is_valid(vec3<int> inds) const;

  inline bool is_factored(int lin) const {
    return _lin2fac.find(lin) != _lin2fac.end();
  }

  double get_h() const { return _h; }

  void visit_neighbors(int lin_center);
  virtual void update_impl(int lin, int * nb, int parent, double & U) = 0;

  struct proxy
  {
    using value_t = double;
    
    proxy(marcher_3d * m): _m {m} {}

    inline value_t get_value(int lin) const {
      return _m->_U[lin];
    }

    inline int get_heap_pos(int lin) const {
      return _m->_heap_pos[lin];
    }

    inline void set_heap_pos(int lin, int pos) {
      _m->_heap_pos[lin] = pos;
    }

  OLIM_PRIVATE:
    marcher_3d * _m {nullptr};
  };

  int _size;

  heap<int, proxy> _heap;

  double * _U {nullptr};
  double const * _s {nullptr};
  state * _state {nullptr};
  int * _heap_pos {nullptr};
  double _h {1};
  vec3<int> _dims;

  // TODO: this is a quick hack just to get this working for the time
  // being.
  std::unordered_map<int, fac_src_3d const *> _lin2fac;
};

#include "marcher_3d.impl.hpp"
