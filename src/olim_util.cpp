#include "olim_util.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

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

template <>
double rhr<2>(double const * p0, double const * dp, double u0, double u1,
              double s_est, double h) {
  check_params(u0, u1, h);
  assert(s_est >= 0);

  double alpha = (u1 - u0)/s_est, alpha_sq = alpha*alpha,
    p0normsq = p0[0]*p0[0] + p0[1]*p0[1],
    d0dotdp = p0[0]*dp[0] + p0[1]*dp[1],
    dpnormsq = dp[0]*dp[0] + dp[1]*dp[1],
    a = (alpha_sq - 1)*p0normsq,
    b = 2*p0dotdp*(alpha_sq - dpnormsq),
    c = alpha_sq*p0normsq - p0dotdp*p0dotdp,
    disc = b*b - 4*a*c;

  assert(disc >= 0);

  double lam = (-b + sqrt(disc))/(2*a),
    plam[2] = {p0[0] + lam*dp[0], p0[1] + lam*dp[1]},
    plamnorm = sqrt(plam[0]*plam[0] + plam[1]*plam[1]),
    plamdotdp = plam[0]*dp[0] + plam[1]*dp[1];

  if (fabs(alpha*plamnorm + plamdotdp) < 1e-13) {
    return (1 - lam)*u0 + lam*u1 + s_est*plamnorm;
  }

  lam = -(b - sqrt(disc))/(2*a);
  plam[0] = p0[0] + lam*dp[0];
  plam[1] = p0[1] + lam*dp[1];
  plamnorm = sqrt(plam[0]*plam[0] + plam[1]*plam[1]);
  plamdotdp = plam[0]*dp[0] + plam[1]*dp[1];

  assert(fabs(alpha*plamnorm + plamdotdp) < 1e-13);

  return (1 - lam)*u0 + lam*u1 + s_est*plamnorm;
}

template <>
double rhr<3>(double const * p0, double const * dp, double u0, double u1,
              double s_est, double h) {
  check_params(u0, u1, h);
  assert(s_est >= 0);

  double alpha = (u1 - u0)/s_est, alpha_sq = alpha*alpha,
    p0normsq = p0[0]*p0[0] + p0[1]*p0[1] + p0[2]*p0[2],
    p0dotdp = p0[0]*dp[0] + p0[1]*dp[1] + p0[2]*dp[2],
    dpnormsq = dp[0]*dp[0] + dp[1]*dp[1] + dp[2]*dp[2],
    a = (alpha_sq - 1)*p0normsq, b = 2*p0dotdp*(alpha_sq - dpnormsq),
    c = alpha_sq*p0normsq - p0dotdp*p0dotdp, disc = b*b - 4*a*c;

  assert(disc >= 0);

  double lam = (-b + sqrt(disc))/(2*a),
    plam[3] = {p0[0] + lam*dp[0], p0[1] + lam*dp[1], p0[2] + lam*dp[2]},
    plamnorm = sqrt(plam[0]*plam[0] + plam[1]*plam[1] + plam[2]*plam[2]),
    plamdotdp = plam[0]*dp[0] + plam[1]*dp[1] + plam[2]*dp[2];

  if (fabs(alpha*plamnorm + plamdotdp) < 1e-13) {
    return (1 - lam)*u0 + lam*u1 + s_est*plamnorm;
  }

  lam = -(b - sqrt(disc))/(2*a);
  plam[0] = p0[0] + lam*dp[0];
  plam[1] = p0[1] + lam*dp[1];
  plam[2] = p0[2] + lam*dp[2];
  plamnorm = sqrt(plam[0]*plam[0] + plam[1]*plam[1] + plam[2]*plam[2]);
  plamdotdp = plam[0]*dp[0] + plam[1]*dp[1] + plam[2]*dp[2];

  assert(fabs(alpha*plamnorm + plamdotdp) < 1e-13);

  return (1 - lam)*u0 + lam*u1 + s_est*plamnorm;
}

/**
 * Adjacent triangle update with constant quadrature rule (used by
 * olim8_rhr and olim8_mp0c).
 */
double rhr_adj(double u0, double u1, double s_est, double h, double * lam) {
  check_params(u0, u1, h);
  assert(s_est >= 0);

  double c = (u0 - u1)/(s_est*h);
  // TODO: use copysign for next line instead?
  double _lam = 0.5 + (c > 0 ? 1 : -1)*std::fabs(c)/(2*sqrt(2 - c*c));
  if (lam != nullptr) {
    *lam = _lam;
  }
  return (1 - _lam)*u0 + _lam*u1 + s_est*h*sqrt(2*_lam*(_lam - 1) + 1);
}

/**
 * Diagonal triangle update with constant quadrature rule (used by
 * olim8_rhr and olim8_mp0c).
 */
double rhr_diag(double u0, double u1, double s_est, double h) {
  check_params(u0, u1, h);
  assert(s_est >= 0);
  
  // TODO: make this one as simple as the adj rule
  double c = std::fabs(u0 - u1)/(s_est*h);
  double sgn = u0 > u1 ? 1 : -1;
  double lam = std::max(0.0, std::min(1.0, sgn*c/sqrt(1 - c*c)));
  return (1 - lam)*u0 + lam*u1 + s_est*h*sqrt(lam*lam + 1);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
