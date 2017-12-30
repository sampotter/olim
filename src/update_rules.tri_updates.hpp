#ifndef __UPDATE_RULES_TRI_UPDATES_HPP__
#define __UPDATE_RULES_TRI_UPDATES_HPP__

#include <type_traits>

#include "speed_est.hpp"
#include "update_rules.utils.hpp"

namespace update_rules {
  template <class speed_est, char degree>
  struct tri_updates: public speed_est {
    using speed_est::speed_est;

    template <char p0, char p1>
    double tri(
      double u0, double u1, double s, double s0, double s1, double h,
      ffvec<p0>, ffvec<p1>) const;

  private:
    template <char p0, char p1>
    double tri_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      ffvec<p0>, ffvec<p1>, std::integral_constant<char, 0>) const;
    template <char p0, char p1>
    double tri_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      ffvec<p0>, ffvec<p1>, std::integral_constant<char, 1>) const;
  };

  using mp0_tri_updates = tri_updates<mp_speed_est, 0>;
  using mp1_tri_updates = tri_updates<mp_speed_est, 1>;
  using rhr_tri_updates = tri_updates<rhr_speed_est, 0>;
}

#include "update_rules.tri_updates.impl.hpp"

#endif // __UPDATE_RULES_TRI_UPDATES_HPP__
