#pragma once

#include <src/config.hpp>

// TODO: remove these
#include <functional>
#include <unordered_map>

#include "heap.hpp"
#include "slow.hpp"
#include "typedefs.h"

struct fac_src
{
  fac_src(double i, double j, double s): i {i}, j {j}, s {s} {}

  double i, j, s;
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

  marcher(int height, int width, double h, no_slow_t const &);
  marcher(int height, int width, double h, double const * s_cache);
  marcher(int height, int width, double h = 1,
          std::function<double(double, double)> s = static_cast<slow2>(s0),
          double x0 = 0.0, double y0 = 0.0);
  virtual ~marcher();

  void init();

  void run();

  void add_boundary_node(int i, int j, double value = 0.0);
  void add_boundary_nodes(int const * i, int const * j, double const * U, int num);
  void add_boundary_nodes(std::tuple<int, int, double> const * nodes, int num);

  void add_boundary_node(double x, double y, double s, double value = 0.0);

  void set_fac_src(int i, int j, fac_src const * src);

  double get_s(int i, int j) const;
  double get_value(int i, int j) const;

  int get_height() const { return _height; }
  int get_width() const { return _width; }

  inline double const * get_s_cache_data() const {
    return _s_cache;
  }

OLIM_PROTECTED:
  inline int linear_index(int i, int j) const {
    return _width*i + j;
  }

  inline int get_i(int lin) const {
    return lin/_width;
  }

  inline int get_j(int lin) const {
    return lin % _width;
  }

  bool in_bounds(int i, int j) const;
  bool in_bounds(double i, double j) const;
  bool is_valid(int i, int j) const;

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
  int _height;
  int _width;

  // TODO: this is a quick hack just to get this working for the time
  // being.
  std::unordered_map<int, fac_src const *> _lin2fac;
};

#include "marcher.impl.hpp"
