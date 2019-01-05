#ifndef __UPDATES_TETRA_HPP__
#define __UPDATES_TETRA_HPP__

#include "cost_funcs.hpp"
#include "updates.common.hpp"
#include "updates.utils.hpp"

namespace updates {

// TODO: maybe we can remove the struct and replace it with a function
// and its specializations?

// TODO: can probably also remove the cost_functor thing
// entirely... Not sure about whether this is worth it or not

template <cost_func F, int n>
struct tetra
{
  void operator()(
    double const * p0, double const * p1, double const * p2,
    double u0, double u1, double u2, double s,
    double s0, double s1, double s2, double h,
    info<2> & info) const;

  void operator()(
    double const * p0, double const * p1, double const * p2,
    double u0, double u1, double u2, double s,
    double s0, double s1, double s2, double h,
    double const * p_fac, double s_fac,
    info<2> & info) const;
};

template <cost_func F, int n, int p0, int p1, int p2>
struct tetra_bv
{
  void operator()(
    double u0, double u1, double u2, double s,
    double s0, double s1, double s2, double h,
    info<2> & info) const;

  // TODO: we aren't actually using this overload yet...
  // template <int n, int p0, int p1, int p2>
  // update_info<2> operator()(F_fac_wkspc<F, 2> & fw) const;
};

}

#include "updates.tetra.impl.hpp"

#endif // __UPDATES_TETRA_HPP__
