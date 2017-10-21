#ifndef __OLIM_RECT_UPDATE_RULES_HPP__
#define __OLIM_RECT_UPDATE_RULES_HPP__

#include "update_rules.line_updates.hpp"

struct olim_rect_update_rules: public update_rules::line_updates {
  double tri11(double u0, double u1, double s, double h) const;
  double tri12(double u0, double u1, double s, double h) const;
  double tri13(double u0, double u1, double s, double h) const;
  double tri22(double u0, double u1, double s, double h) const;
  double tri23(double u0, double u1, double s, double h) const;
  double tetra111(double u0, double u1, double u2, double s, double h) const;
  double tetra122(double u0, double u1, double u2, double s, double h) const;
  double tetra123(double u0, double u1, double u2, double s, double h) const;
  double tetra222(double u0, double u1, double u2, double s, double h) const;
};

#endif // __OLIM_RECT_UPDATE_RULES_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
