#include "olim8_rhr.hpp"

#include <cmath>

#include "olim_util.hpp"

double olim8_rhr_update_rules::adj1pt(double u0, double s, double s0,
                                      double h) const {
  (void) s0;
  return u0 + s*h;
}

double olim8_rhr_update_rules::adj2pt(double u0, double u1, double s, double s0,
                                      double s1, double h) const {
  (void) s0;
  (void) s1;
  return rhr_adj(u0, u1, s, h);
}

double olim8_rhr_update_rules::diag1pt(double u0, double s, double s0,
                                       double h) const {
  (void) s0;
  return u0 + s*h*std::sqrt(2);
}

double olim8_rhr_update_rules::diag2pt(double u0, double u1, double s,
                                       double s0, double s1, double h) const {
  (void) s0;
  (void) s1;
  return rhr_diag(u0, u1, s, h);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
