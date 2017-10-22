#include "olim8_rhr.hpp"

#include <cmath>

#include "common.defs.hpp"
#include "olim_util.hpp"

double olim8_rhr_update_rules::tri11(double u0, double u1, double s, double s0,
                                     double s1, double h) const {
  (void) s0;
  (void) s1;
  return rhr_adj(u0, u1, s, h);
}

double olim8_rhr_update_rules::tri12(double u0, double u1, double s,
                                     double s0, double s1, double h) const {
  (void) s0;
  (void) s1;
  return rhr_diag(u0, u1, s, h);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
