#pragma once

#include <src/config.hpp>

// TODO: remove these
#include <unordered_map>

#include "../base_marcher.hpp"
#include "../fac.hpp"
#include "../state.hpp"
#include "../vec.hpp"

#include "slow.hpp"

namespace eikonal {

template <
  class base,
  int n,
  int num_nb,
  ordering ord
>
struct marcher: public base_marcher<marcher<base, n, num_nb, ord>, n>
{
  using base_marcher_t = base_marcher<marcher<base, n, num_nb, ord>, n>;

  using fac_src_t = fac_src<n>;

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

  void add_bd(int const * inds);
  void add_bd(ivec inds);

  void add_free(int const * inds);
  void add_free(ivec inds);

  void set_fac_src(int * inds, fac_src<n> const * src);
  void set_fac_src(ivec inds, fac_src<n> const * src);

  double get_s(ivec inds) const;

  inline double * get_s_ptr() const {
    return _s;
  }

  inline void set_s_ptr(double * s) {
    _s = s;
  }

OLIM_PROTECTED:
  inline bool is_factored(int lin) const {
    return _lin2fac.find(lin) != _lin2fac.end();
  }

  double get_h() const { return _h; }

  void visit_neighbors(int lin);

  double * _s {nullptr};
  double _h {1};

  int _child_offset[num_nb][num_nb];

  // TODO: a quick hack just to get this working for now
  std::unordered_map<int, fac_src<n> const *> _lin2fac;
};

}

#include "marcher.impl.hpp"
