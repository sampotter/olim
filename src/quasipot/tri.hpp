#pragma once

#include "../update_info.hpp"
#include "../vec.hpp"

namespace quasipot {

// TODO: an idea for how to reduce a lot of boilerplate here:
// - replace (U0, U1), (b0, b1), and (x0, x1) with vectors/matrices
// - add a template parameter `d`
// - most everything else below follows...
// - except, we will need to be a little careful of how we
//   implement the solver...
// - can we solve explicitly without using a rootfinder?
// - it may be hard to beat having special overloads for one-point and
//   two-point updates

template <int n>
struct tri
{
  using vec_t = vec<double, n>;

  void operator()(
    vec_t const & x0, vec_t const & x1,
    double U0, double U1,
    vec_t const & b0, vec_t const & b1,
    update_info<1> & info
    ) const
  {
    double dU = U1 - U0, xt_norm, bt_norm;
    vec_t db = b1 - b0, dx = x1 - x0, xt, bt;

    auto const dF = [&] (double t) {
      xt = x0 + t*dx;
      bt = b0 + t*db;
      xt_norm = xt.norm2();
      bt_norm = bt.norm2();
      return dU + (xt_norm*db - bt_norm*dx)*(bt/bt_norm - xt/xt_norm);
    };

    double t;
    hybrid_status status;
    std::tie(t, status) = hybrid(dF, 0., 1.);

    if (status == hybrid_status::DEGENERATE) {
      double F0 = U0 + b0.norm2()*x0.norm2() - b0*x0;
      double F1 = U1 + b1.norm2()*x1.norm2() - b1*x1;
      info.lambda[0] = F0 < F1 ? 0 : 1;
      info.value = fmin(F0, F1);
    } else {
      info.lambda[0] = t;
      // TODO: try this but without recomputing---we might be okay
      // just using the last computations from dF
      xt = x0 + t*dx;
      bt = b0 + t*db;
      info.value = U0 + t*dU + bt.norm2()*xt.norm2() - bt*xt;
    }
  }
};

}
