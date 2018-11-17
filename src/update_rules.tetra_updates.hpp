#ifndef __UPDATE_RULES_TETRA_UPDATES_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_HPP__

#include "cost_funcs.hpp"
#include "update_rules.utils.hpp"

namespace update_rules {

template <cost_func F_>
struct tetra
{
  update_info<2> operator()(
    double const * p0, double const * p1, double const * p2,
    double u0, double u1, double u2, double s,
    double s0, double s1, double s2, double h) const;

  update_info<2> operator()(
    double const * p0, double const * p1, double const * p2,
    double u0, double u1, double u2, double s,
    double s0, double s1, double s2, double h,
    double const * p_fac, double s_fac) const;
};

using tetra_mp0 = tetra<cost_func::mp0>;
using tetra_mp1 = tetra<cost_func::mp1>;
using tetra_rhr = tetra<cost_func::rhr>;

template <cost_func F, char p0, char p1, char p2>
struct tetra_bv {
  using wkspc = F_wkspc<F, 2>;
  using fac_wkspc = F_fac_wkspc<F, 2>;
  update_info<2> operator()(wkspc const & w) const;
  update_info<2> operator()(fac_wkspc const & fw) const;
};

template <char p0, char p1, char p2>
using tetra_bv_mp0 = tetra_bv<cost_func::mp0, p0, p1, p2>;

template <char p0, char p1, char p2>
using tetra_bv_mp1 = tetra_bv<cost_func::mp1, p0, p1, p2>;

template <char p0, char p1, char p2>
using tetra_bv_rhr = tetra_bv<cost_func::rhr, p0, p1, p2>;

}

#include "update_rules.tetra_updates.impl.hpp"

#endif // __UPDATE_RULES_TETRA_UPDATES_HPP__
