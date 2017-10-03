#include "olim_util.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>

#ifdef EIKONAL_DEBUG
void check_params(double u0, double u1, double h, double s) {
  assert(u0 >= 0);
  assert(u1 >= 0);
  assert(!std::isinf(u0));
  assert(!std::isinf(u1));
  assert(h > 0);
  assert(s >= 0);
}
#endif

#ifdef EIKONAL_DEBUG
void check_params(double u0, double u1, double u2, double h, double s) {
  assert(u0 >= 0);
  assert(u1 >= 0);
  assert(u2 >= 0);
  assert(!std::isinf(u0));
  assert(!std::isinf(u1));
  assert(!std::isinf(u2));
  assert(h > 0);
  assert(s >= 0);
}
#endif

/**
 * Adjacent triangle update with constant quadrature rule (used by
 * olim8_rhr and olim8_mp0c).
 */
double rhr_adj(double u0, double u1, double s, double h, double * lam) {
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif

  double c = (u0 - u1)/(s*h);
  // TODO: use copysign for next line instead?
  double _lam = 0.5 + (c > 0 ? 1 : -1)*std::fabs(c)/(2*sqrt(2 - c*c));
  if (lam != nullptr) {
    *lam = _lam;
  }
  return (1 - _lam)*u0 + _lam*u1 + s*h*sqrt(2*_lam*(_lam - 1) + 1);
}

/**
 * Diagonal triangle update with constant quadrature rule (used by
 * olim8_rhr and olim8_mp0c).
 */
double rhr_diag(double u0, double u1, double s, double h) {
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  
  // TODO: make this one as simple as the adj rule
  double sh = s*h, c = std::fabs(u0 - u1)/sh;
  double sgn = u0 > u1 ? 1 : -1;
  double lam = std::max(0.0, std::min(1.0, sgn*c/sqrt(1 - c*c)));
  return (1 - lam)*u0 + lam*u1 + sh*sqrt(lam*lam + 1);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
