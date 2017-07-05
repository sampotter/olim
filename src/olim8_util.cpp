#include "olim8_util.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

#include <gsl/gsl_poly.h>

static void check_params(double u0, double u1, double h) {
  (void) u0;
  (void) u1;
  (void) h;
  assert(h > 0);
  assert(u0 >= 0);
  assert(u1 >= 0);
  assert(!std::isinf(u0));
  assert(!std::isinf(u1));
}

/**
 * Adjacent triangle update with constant quadrature rule (used by
 * olim8_rhr and olim8_mp0).
 */
double rhr_adj(double u0, double u1, double s_est, double h) {
  check_params(u0, u1, h);
  assert(s_est >= 0);

  double c = (u0 - u1)/(s_est*h);
  // TODO: use copysign for next line instead?
  double lam = 0.5 + (c > 0 ? 1 : -1)*std::fabs(c)/(2*std::sqrt(2 - c*c));
  return (1 - lam)*u0 + lam*u1 + s_est*h*sqrt(2*lam*(lam - 1) + 1);
}

/**
 * Diagonal triangle update with constant quadrature rule (used by
 * olim8_rhr and olim8_mp0).
 */
double rhr_diag(double u0, double u1, double s_est, double h) {
  check_params(u0, u1, h);
  assert(s_est >= 0);
  
  // TODO: make this one as simple as the adj rule
  double c = std::fabs(u0 - u1)/(s_est*h);
  double sgn = u0 > u1 ? 1 : -1;
  double lam = std::max(0.0, std::min(1.0, sgn*c/std::sqrt(1 - c*c)));
  return (1 - lam)*u0 + lam*u1 + s_est*h*sqrt(lam*lam + 1);
}

double mp1_solve(double * a, double u0, double u1, double s0, double s1,
                 double h) {
  /**
   * Solve quartic
   */
  double z[8];
  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(5);
  gsl_poly_complex_solve(a, 5, w, z);
  gsl_poly_complex_workspace_free(w);

  /**
   * Extract a reasonable (i.e. nonnegative) value from the real
   * roots---return inf otherwise.
   */
  double lam = -1, uhat = std::numeric_limits<double>::infinity();
  for (int i = 0; i < 4; ++i) {
    if (z[2*i + 1] != 0) {
      continue;
    }
    lam = z[2*i];
    if (0 <= lam && lam <= 1) {
      double val = (1 - lam)*u0 + lam*u1 +
        h*((1 - lam)*s0 + lam*s1)*std::sqrt(2*lam*lam - 2*lam + 1);
      if (val >= 0) {
        uhat = std::min(uhat, val);
      }
    }
  }
  return uhat;
}

double mp1_adj(double u0, double u1, double s0, double s1, double h) {
  check_params(u0, u1, h);
  assert(s0 >= 0);
  assert(s1 >= 0);

  double alpha_sq = std::pow((u0 - u1)/h, 2);
  double s0_sq = s0*s0;
  double ds = s1 - s0;
  double ds_sq = ds*ds;
  double a[] = {
    s0_sq - alpha_sq - 2*s0*ds + ds_sq,
    -4*s0_sq + 2*alpha_sq + 10*s0*ds - 6*ds_sq,
    4*s0_sq - 2*alpha_sq - 20*s0*ds + 17*ds_sq,
    16*s0*ds - 24*ds_sq,
    16*ds_sq
  };

  return mp1_solve(a, u0, u1, s0, s1, h);
}

double mp1_diag(double u0, double u1, double s0, double s1, double h) {
  check_params(u0, u1, h);
  assert(s0 >= 0);
  assert(s1 >= 0);

  double alpha_sq = std::pow((u0 - u1)/h, 2);
  double ds = s1 - s0;
  double ds_sq = ds*ds;
  double a[] = {ds_sq - alpha_sq, 2*s0*ds, s0*s0 + 4*ds_sq, 4*s0*ds, 4*ds_sq};

  return mp1_solve(a, u0, u1, s0, s1, h);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
