#pragma once

#include <src/config.hpp>

// TODO: remove these
#include <unordered_map>

#include "../base_marcher.hpp"
#include "../fac.h"
#include "../state.hpp"
#include "../vec.hpp"

#include "slow.hpp"

namespace eikonal {

template <
  class derived,
  int n,
  int num_nb,
  ordering ord
>
struct marcher: public base_marcher<marcher<derived, n, num_nb, ord>, n>
{
  using base_marcher_t = base_marcher<marcher<derived, n, num_nb, ord>, n>;

  using fvec = vec<double, n>;
  using ivec = vec<int, n>;
  using uvec = vec<unsigned, n>;

  static constexpr int ndim = n; // TODO: -> get_n()

  static constexpr int get_num_nb() { return num_nb; }
  static constexpr ordering get_ord() { return ord; }

  marcher(ivec dims, double h, no_slow_t const &);
  marcher(ivec dims, double h);
  marcher(ivec dims, double h, double const * s);
  ~marcher();

  using base_marcher_t::add_src;
  void add_src(fvec coords, double s, double U = 0.0);

  void factor(int * inds, fac_src const * src);
  void factor(ivec inds, fac_src const * src);

  double get_s(ivec inds) const;

  inline double * get_s_ptr() const {
    return _s;
  }

  inline void set_s_ptr(double * s) {
    _s = s;
  }

  inline bool is_factored(int lin) const {
    return _fac_srcs.find(lin) != _fac_srcs.end();
  }

  double get_h() const { return _h; }

  void visit_neighbors(int lin);

  double * _s {nullptr};
  double _h {1};

  int _child_offset[num_nb][num_nb];

  // TODO: a quick hack just to get this working for now
  std::unordered_map<int, fac_src const *> _fac_srcs;
};

}

#include "marcher.impl.hpp"
