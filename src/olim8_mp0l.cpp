#include "olim8_mp0l.hpp"

#include <cmath>

#include "olim_util.hpp"

double olim8_mp0l_update_rules::adj1pt(double u0, double s, double s0,
                                       double h) const {
  return u0 + h*(s + s0)/2;
}

double olim8_mp0l_update_rules::adj2pt(double u0, double u1, double s,
                                       double s0, double s1, double h) const {
  return mp0l_adj(u0, u1, s, s0, s1, h);
}

double olim8_mp0l_update_rules::diag1pt(double u0, double s, double s0,
                                        double h) const {
  return u0 + h*(s + s0)*std::sqrt(2)/2;
}

double olim8_mp0l_update_rules::diag2pt(double u0, double u1, double s,
                                        double s0, double s1, double h) const {
  return mp0l_diag(u0, u1, s, s0, s1, h);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
