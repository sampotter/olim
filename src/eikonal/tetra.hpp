#pragma once

#include "func.hpp"

#include "../update_info.hpp"
#include "../vec.hpp"

namespace eikonal {

// TODO: maybe we can remove the struct and replace it with a function
// and its specializations?
//
// TODO: there are no more specializations---we can definitely do
// this. Also remove the impl file and just move everything in here at
// the same time.

// TODO: can probably also remove the cost_functor thing
// entirely... Not sure about whether this is worth it or not

template <cost_func F, int n>
struct tetra {
  void operator()(cost_functor<F, n, 2> & func, update_info<2> & info) const;
  void operator()(cost_functor_fac<F, n, 2> & func, update_info<2> & info) const;
};

template <cost_func F, int n, int p0, int p1, int p2>
struct tetra_bv
{
  void operator()(cost_functor_bv<F, n, p0, p1, p2> & func, update_info<2> & info) const;

  // TODO: we aren't actually using this overload yet...
  // template <int n, int p0, int p1, int p2>
  // update_info<2> operator()(F_fac_wkspc<F, 2> & fw) const;
};

}

#include "tetra.impl.hpp"
