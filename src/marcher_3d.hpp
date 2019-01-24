#ifndef __MARCHER_3D_HPP__
#define __MARCHER_3D_HPP__

// TODO: try to remove this
#include <functional>
#include <unordered_map>

#include "heap.hpp"
#include "speed_funcs.hpp"
#include "typedefs.h"

struct fac_src_3d
{
  fac_src_3d(double i, double j, double k, double s):
    i {i}, j {j}, k {k}, s {s} {}

  double i, j, k, s;
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
  marcher_3d(int height, int width, int depth, double h,
             no_speed_func_t const &);
  marcher_3d(int height, int width, int depth, double h = 1,
             std::function<double(double, double, double)> speed =
             static_cast<speed_func_3d>(default_speed_func),
             double x0 = 0.0, double y0 = 0.0, double z0 = 0.0);
  marcher_3d(int height, int width, int depth, double h,
             double const * s_cache);
  virtual ~marcher_3d();

  void init();

  void run();

  void add_boundary_node(int i, int j, int k, double value = 0.0);
  void add_boundary_nodes(int const * i, int const * j, int const * k,
                          double const * U, int num);

  void add_boundary_node(
    double x, double y, double z, double s, double value = 0.0);

  void set_fac_src(int i, int j, int k, fac_src_3d const * src);

  double get_speed(int i, int j, int k) const;
  double get_value(int i, int j, int k) const;

  int get_height() const { return _height; }
  int get_width() const { return _width; }
  int get_depth() const { return _depth; }

  inline double const * get_s_cache_data() const {
    return _s_cache;
  }

OLIM_PROTECTED:
  inline int linear_index(int i, int j, int k) const {
    return _height*(_width*k + j) + i; // column-major
  }

  inline int get_i(int lin) const {
    return lin % _height;
  }

  inline int get_j(int lin) const {
    return lin/_height % _width;
  }

  inline int get_k(int lin) const {
    return lin/(_height*_width);
  }

  bool in_bounds(int i, int j, int k) const;
  bool is_valid(int i, int j, int k) const;

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
  double const * _s_cache {nullptr};
  state * _state {nullptr};
  int * _heap_pos {nullptr};
  double _h {1};
  int _height;
  int _width;
  int _depth;

  // TODO: this is a quick hack just to get this working for the time
  // being.
  std::unordered_map<int, fac_src_3d const *> _lin2fac;
};

#include "marcher_3d.impl.hpp"

#endif // __MARCHER_3D_HPP__
