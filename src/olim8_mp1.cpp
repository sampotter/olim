#include "olim8_mp1.hpp"

#include <cmath>

#include "olim_util.hpp"

double olim8_mp1_update_rules::adj1pt(double u0, double s, double s0,
                                      double h) const {
  return u0 + h*(s + s0)/2;
}

double olim8_mp1_update_rules::adj2pt(double u0, double u1, double s,
                                      double s0, double s1, double h) const {
  double sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
  return sbar0 == sbar1 ? rhr_adj(u0, u1, sbar0, h) :
    mp1_adj(u0, u1, sbar0, sbar1, h);
}

double olim8_mp1_update_rules::diag1pt(double u0, double s, double s0,
                                       double h) const {
  return u0 + h*(s + s0)*std::sqrt(2)/2;
}

double olim8_mp1_update_rules::diag2pt(double u0, double u1, double s,
                                       double s0, double s1, double h) const {
  double sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
  return sbar0 == sbar1 ? rhr_diag(u0, u1, sbar0, h) :
    mp1_diag(u0, u1, sbar0, sbar1, h);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
