#include "olim8_mp0.hpp"

#include <cmath>

#include "common.defs.hpp"
#include "olim_util.hpp"

double olim8_mp0_update_rules::tri11(double u0, double u1, double s,
                                     double s0, double s1, double h) const {
  return rhr_adj(u0, u1, (s + (s0 + s1)/2)/2, h);
}

double olim8_mp0_update_rules::tri12(double u0, double u1, double s,
                                     double s0, double s1, double h) const {
  return rhr_diag(u0, u1, (s + (s0 + s1)/2)/2, h);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
