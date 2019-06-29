#pragma once

#include "../update_info.hpp"
#include "../vec.hpp"

namespace quasipot {

template <int n>
struct tetra
{
  using vec_t = vec<double, n>;

  void operator()(
    vec_t const & x0, vec_t const & x1, vec_t const & x2,
    double U0, double U1, double U2,
    vec_t const & b0, vec_t const & b1, vec_t const & b2
    ) const
  {
    mat<double, n, 2> dX {x1 - x0, x2 - x0};
    vec2<double> dU {U1 - U0, U2 - U0};
    mat<double, n, 2> dB {b1 - b0, b2 - b0};

    auto const dF = [&] (vec2<double> const & t) {
      xt = x0 + dU*t;
      bt = b0 + dB*t;
      xt_norm = xt.norm();
      bt_norm = bt.norm();
      return dU + (xt_norm*dB - bt_norm*dX).t()*(bt/bt_norm - xt/xt_norm);
    };

    auto const d2F = [&] (vec2<double> const & t) {
      xt = x0 + dU*t;
      bt = b0 + dB*t;

    };
  }
};

}
