#ifndef __UPDATE_RULES_TETRA_UPDATES_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_HPP__

#include "speed_estimators.hpp"

namespace update_rules {
  template <class speed_estimator>
  struct rect_tetra_updates: public speed_estimator {
    double tetra111(double u0, double u1, double u2, double s,
                    double s0, double s1, double s2, double h) const;
    double tetra122(double u0, double u1, double u2, double s,
                    double s0, double s1, double s2, double h) const;
    double tetra123(double u0, double u1, double u2, double s,
                    double s0, double s1, double s2, double h) const;
    double tetra222(double u0, double u1, double u2, double s,
                    double s0, double s1, double s2, double h) const;
  };

  using rhr_tetra_updates = rect_tetra_updates<rhr_speed_estimator>;
  using mp0_tetra_updates = rect_tetra_updates<mp_speed_estimator>;

  struct mp1_tetra_updates {
    double tetra111(double u0, double u1, double u2, double s,
                    double s0, double s1, double s2, double h) const;
    double tetra122(double u0, double u1, double u2, double s,
                    double s0, double s1, double s2, double h) const;
    double tetra123(double u0, double u1, double u2, double s,
                    double s0, double s1, double s2, double h) const;
    double tetra222(double u0, double u1, double u2, double s,
                    double s0, double s1, double s2, double h) const;
  };
}

#include "update_rules.tetra_updates.impl.hpp"

#endif // __UPDATE_RULES_TETRA_UPDATES_HPP__
