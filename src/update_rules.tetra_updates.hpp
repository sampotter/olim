#ifndef __UPDATE_RULES_TETRA_UPDATES_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_HPP__

#include "speed_est.hpp"
#include "update_rules.utils.hpp"

namespace update_rules {
  template <class speed_est, char degree>
  struct tetra_updates: public speed_est {
    using speed_est::speed_est;

    template <char p0, char p1, char p2>
    double tetra(
      double u0, double u1, double u2, double s,
      double s0, double s1, double s2, double h,
      ffvec<p0>, ffvec<p1>, ffvec<p2>) const;

  private:
    template <char p0, char p1, char p2>
    double tetra_impl(
      double u0, double u1, double u2, double s,
      double s0, double s1, double s2, double h,
      ffvec<p0>, ffvec<p1>, ffvec<p2>, std::integral_constant<char, 0>) const;
    template <char p0, char p1, char p2>
    double tetra_impl(
      double u0, double u1, double u2, double s,
      double s0, double s1, double s2, double h,
      ffvec<p0>, ffvec<p1>, ffvec<p2>, std::integral_constant<char, 1>) const;
  };

  using mp0_tetra_updates = tetra_updates<mp_speed_est, 0>;
  using mp1_tetra_updates = tetra_updates<mp_speed_est, 1>;
  using rhr_tetra_updates = tetra_updates<rhr_speed_est, 0>;
}

#include "update_rules.tetra_updates.impl.hpp"

#endif // __UPDATE_RULES_TETRA_UPDATES_HPP__
