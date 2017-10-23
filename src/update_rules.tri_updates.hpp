#ifndef __UPDATE_RULES_TRI_UPDATES_HPP__
#define __UPDATE_RULES_TRI_UPDATES_HPP__

#include <type_traits>

#include "speed_estimators.hpp"

// TODO: consider cleaning up the template mess below using a macro
//
// TODO: replace the tri functions with a tri<D1, D2> function (one
// benefit of this is that the template parameter is available for use
// in the implementation of these functions, which should help to keep
// bugs at bay)

namespace update_rules {
  template <class speed_estimator, bool is_constrained>
  struct rect_tri_updates: public speed_estimator {
    double tri11(double u0, double u1, double s, double s0, double s1, double h) const;
    double tri12(double u0, double u1, double s, double s0, double s1, double h) const;
    double tri13(double u0, double u1, double s, double s0, double s1, double h) const;
    double tri22(double u0, double u1, double s, double s0, double s1, double h) const;
    double tri23(double u0, double u1, double s, double s0, double s1, double h) const;

  private:
    double
    tri11_impl(double u0, double u1, double s, double s0, double s1, double h,
               std::true_type &&) const;

    double
    tri11_impl(double u0, double u1, double s, double s0, double s1, double h,
               std::false_type &&) const;

    double
    tri12_impl(double u0, double u1, double s, double s0, double s1, double h,
               std::true_type &&) const;

    double
    tri12_impl(double u0, double u1, double s, double s0, double s1, double h,
               std::false_type &&) const;

    double
    tri13_impl(double u0, double u1, double s, double s0, double s1, double h,
               std::false_type &&)
      const;

    double
    tri22_impl(double u0, double u1, double s, double s0, double s1, double h,
               std::false_type &&)
      const;

    double
    tri23_impl(double u0, double u1, double s, double s0, double s1, double h,
               std::false_type &&)
      const;
  };

  template <bool is_constrained>
  using rhr_tri_updates = rect_tri_updates<rhr_speed_estimator, is_constrained>;

  template <bool is_constrained>
  using mp0_tri_updates = rect_tri_updates<mp_speed_estimator, is_constrained>;

  template <bool is_constrained>
  struct mp1_tri_updates {
    double tri11(double u0, double u1, double s, double s0, double s1, double h)
      const;
    double tri12(double u0, double u1, double s, double s0, double s1, double h)
      const;
    double tri13(double u0, double u1, double s, double s0, double s1, double h)
      const;
    double tri22(double u0, double u1, double s, double s0, double s1, double h)
      const;
    double tri23(double u0, double u1, double s, double s0, double s1, double h)
      const;

  private:
    double tri11_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      std::true_type &&) const;

    double tri11_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      std::false_type &&) const;

    double tri12_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      std::true_type &&) const;

    double tri12_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      std::false_type &&) const;

    double tri13_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      std::false_type &&) const;

    double tri22_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      std::false_type &&) const;

    double tri23_impl(
      double u0, double u1, double s, double s0, double s1, double h,
      std::false_type &&) const;
  };
}

#include "update_rules.tri_updates.impl.hpp"

#endif // __UPDATE_RULES_TRI_UPDATES_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
