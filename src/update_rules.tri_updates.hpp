#ifndef __UPDATE_RULES_TRI_UPDATES_HPP__
#define __UPDATE_RULES_TRI_UPDATES_HPP__

#include <type_traits>

#include "speed_estimators.hpp"

namespace update_rules {
  template <class speed_estimator, bool is_constrained>
  struct rect_tri_updates: public speed_estimator {
    std::enable_if_t<!is_constrained, double>
    tri11(double u0, double u1, double s, double s0, double s1, double h) const;

    std::enable_if_t<!is_constrained, double>
    tri12(double u0, double u1, double s, double s0, double s1, double h) const;

    std::enable_if_t<!is_constrained, double>
    tri13(double u0, double u1, double s, double s0, double s1, double h) const;

    std::enable_if_t<!is_constrained, double>
    tri22(double u0, double u1, double s, double s0, double s1, double h) const;

    std::enable_if_t<!is_constrained, double>
    tri23(double u0, double u1, double s, double s0, double s1, double h) const;
  };

  template <bool is_constrained>
  using rhr_tri_updates = rect_tri_updates<rhr_speed_estimator, is_constrained>;

  template <bool is_constrained>
  using mp0_tri_updates = rect_tri_updates<mp_speed_estimator, is_constrained>;

  // TODO: struct mp1_tri_updates
}

#include "update_rules.tri_updates.impl.hpp"

#endif // __UPDATE_RULES_TRI_UPDATES_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
