#ifndef __OLIM18_RHR_IMPL_HPP__
#define __OLIM18_RHR_IMPL_HPP__

#include <cmath>
#include <cstdio>

#include "olim_util.hpp"

template <class rootfinder>
double olim18_rhr_update_rules<rootfinder>::line1(
  double u0, double s, double s0, double h) const
{
  (void) s0;
  return u0 + s*h;
}

template <class rootfinder>
double olim18_rhr_update_rules<rootfinder>::line2(
  double u0, double s, double s0, double h) const
{
  (void) s0;
  return u0 + s*h*sqrt(2);
}

template <class rootfinder>
double olim18_rhr_update_rules<rootfinder>::tri12(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
  return rhr_diag(u0, u1, s, h);
}

template <class rootfinder>
double olim18_rhr_update_rules<rootfinder>::tri22(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
  return rhr_3d_22(u0, u1, s, h);
}

template <class rootfinder>
double olim18_rhr_update_rules<rootfinder>::tetra122(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  (void) s0;
  (void) s1;
  (void) s2;

  double sh = s*h, du1 = u1 - u0, du2 = u2 - u0;
  double alpha1 = fabs(du1/sh), alpha1_sq = alpha1*alpha1;
  double alpha2 = fabs(du2/sh), alpha2_sq = alpha2*alpha2;

  assert(fabs(1 - 2*alpha1_sq*alpha2_sq) > 1e-13);
  double denom = sqrt(1 - 2*alpha1_sq*alpha2_sq);

  double T = std::numeric_limits<double>::infinity();
  double lam1, lam2, l;

  lam1 = sqrt(alpha1_sq + (alpha2_sq + alpha1_sq)*alpha2_sq)/denom;
  lam2 = sqrt((alpha1_sq - alpha2_sq)*alpha2_sq + alpha2_sq)/denom;
  l = sqrt(1 + lam1*lam1 + lam2*lam2);
  if (fabs(du1*l/sh + lam1) < 1e-13 && fabs(du2*l/sh + lam2) < 1e-13) {
    T = std::min(T, (1 - lam1 - lam2)*u0 + lam1*u1 + lam2*u2 + sh*l);
  }

  lam1 = -lam1;
  l = sqrt(1 + lam1*lam1 + lam2*lam2);
  if (fabs(du1*l/sh + lam1) < 1e-13 && fabs(du2*l/sh + lam2) < 1e-13) {
    T = std::min(T, (1 - lam1 - lam2)*u0 + lam1*u1 + lam2*u2 + sh*l);
  }

  lam2 = -lam2;
  l = sqrt(1 + lam1*lam1 + lam2*lam2);
  if (fabs(du1*l/sh + lam1) < 1e-13 && fabs(du2*l/sh + lam2) < 1e-13) {
    T = std::min(T, (1 - lam1 - lam2)*u0 + lam1*u1 + lam2*u2 + sh*l);
  }

  lam1 = -lam1;
  l = sqrt(1 + lam1*lam1 + lam2*lam2);
  if (fabs(du1*l/sh + lam1) < 1e-13 && fabs(du2*l/sh + lam2) < 1e-13) {
    T = std::min(T, (1 - lam1 - lam2)*u0 + lam1*u1 + lam2*u2 + sh*l);
  }

  return T;
}

template <class rootfinder>
double olim18_rhr_update_rules<rootfinder>::tetra222(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  (void) s0;
  (void) s1;
  (void) s2;

  double sh_scaled = s*h*2*sqrt(2), du1 = u1 - u0, du2 = u2 - u0;
  double alpha1 = fabs(du1/sh_scaled), alpha2 = fabs(du2/sh_scaled);
  double alpha1_sq = alpha1*alpha1, alpha2_sq = alpha2*alpha2;
  double A1 = 2 - alpha1_sq, B1 = 1 - 2*alpha1_sq, C1 = 1 - alpha1_sq;
  double A2 = 2 - alpha2_sq, B2 = 1 - 2*alpha2_sq, C2 = 1 - alpha2_sq;
  double AQ1[6] = {2*A1, 2*A1, B1, -2*A1, -2*C1, C1};
  double AQ2[6] = {B2, 2*A2, 2*A2, -2*C2, -2*A1, C2};

  int n = 0;
  double isects[8];
  this->intersect_conics(AQ1, AQ2, isects, n);

  double T = std::numeric_limits<double>::infinity();
  double lam0, lam1, lam2, l;
  for (int i = 0; i < n; i += 2) {
    lam1 = isects[i];
    lam2 = isects[i + 1];
    lam0 = 1 - lam1 - lam2;
    l = sqrt(lam0*lam0 + lam1*lam1 + lam2*lam2);
    if (fabs(2*lam1 + lam2 - 1 + du1*l/sh_scaled) < 1e-13 &&
        fabs(2*lam1 + lam2 - 1 + du1*l/sh_scaled) < 1e-13) {
      T = std::min(T, lam0*u0 + lam1*u1 + lam2*u2 + sh_scaled*l);
    }
  }
  return T;
}

#endif // __OLIM18_RHR_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
