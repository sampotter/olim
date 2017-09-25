#ifndef __OLIM_UPDATE_RULES_HPP__
#define __OLIM_UPDATE_RULES_HPP__

template <class rootfinder>
struct olim3d_rhr_update_rules: rootfinder {
  double line1(double u0, double s, double s0, double h) const;
  double line2(double u0, double s, double s0, double h) const;
  double line3(double u0, double s, double s0, double h) const;
  double tri11(double u0, double u1, double s,
               double s0, double s1, double h) const;
  double tri12(double u0, double u1, double s,
               double s0, double s1, double h) const;
  double tri13(double u0, double u1, double s,
               double s0, double s1, double h) const;
  double tri22(double u0, double u1, double s,
               double s0, double s1, double h) const;
  double tri23(double u0, double u1, double s,
               double s0, double s1, double h) const;
  double tetra111(double u0, double u1, double u2, double s,
                  double s0, double s1, double s2, double h) const;
  double tetra122(double u0, double u1, double u2, double s,
                  double s0, double s1, double s2, double h) const;
  double tetra123(double u0, double u1, double u2, double s,
                  double s0, double s1, double s2, double h) const;
  double tetra222(double u0, double u1, double u2, double s,
                  double s0, double s1, double s2, double h) const;
};

#include "olim_update_rules.impl.hpp"

#endif // __OLIM_UPDATE_RULES_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
