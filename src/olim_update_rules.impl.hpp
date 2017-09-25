#ifndef __OLIM_UPDATE_RULES_IMPL_HPP__
#define __OLIM_UPDATE_RULES_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.defs.hpp"
#include "olim_update_rules.hpp"
#include "olim_util.hpp"

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::line1(
  double u0, double s, double s0, double h) const
{
  (void) s0;
#if PRINT_UPDATES
  double tmp = u0 + s*h;
  printf("line1(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h;
#endif
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::line2(
  double u0, double s, double s0, double h) const
{
  (void) s0;
#if PRINT_UPDATES
  double tmp = u0 + s*h*sqrt2;
  printf("line2(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h*sqrt2;
#endif
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::line3(
  double u0, double s, double s0, double h) const
{
  (void) s0;
#if PRINT_UPDATES
  double tmp = u0 + s*h*sqrt3;
  printf("line3(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h*sqrt3;
#endif
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tri11(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
#if PRINT_UPDATES
  double tmp = rhr_adj(u0, u1, s, h);
  printf("tri11(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, tmp);
  return tmp;
#else
  return rhr_adj(u0, u1, s, h);
#endif
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tri12(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
#if PRINT_UPDATES
  double tmp = rhr_diag(u0, u1, s, h);
  printf("tri12(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, tmp);
  return tmp;
#else
  return rhr_diag(u0, u1, s, h);
#endif
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tri13(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, alpha = fabs(u1 - u0)/sh, sgn = u0 > u1 ? 1 : -1;
  assert(2 > alpha*alpha);
  double lam = std::max(0.0, std::min(1.0, sgn*alpha/sqrt(2*(2 - alpha*alpha))));
#if PRINT_UPDATES
  double tmp = (1 - lam)*u0 + lam*u1 + sh*sqrt(1 + 2*lam*lam);
  printf("tri13(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, tmp);
  return tmp;
#else
  return (1 - lam)*u0 + lam*u1 + sh*sqrt(1 + 2*lam*lam);
#endif
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tri22(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double du = u1 - u0, sh = s*h, alpha = du/sh, alpha_sq = alpha*alpha;
  double disc = 1 - 2*(2*alpha_sq - 1)/(alpha_sq - 2);
  assert(disc >= 0);
  double lam = std::max(
    0.0, std::min(1.0, (1 + (u0 > u1 ? 1 : -1)*std::sqrt(disc))/2));
  double llam = sqrt(2*(lam*(lam - 1) + 1));
#if PRINT_UPDATES
  double tmp = (1 - lam)*u0 + lam*u1 + sh*llam;
  printf("tri22(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, tmp);
  return tmp;
#else
  return (1 - lam)*u0 + lam*u1 + sh*llam;
#endif
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tri23(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, alpha = fabs(u1 - u0)/sh, sgn = u0 > u1 ? 1 : -1;
  assert(1 >= alpha*alpha);
  double lam = std::max(
    0.0, std::min(1.0, sgn*sqrt2*alpha/sqrt(1 - alpha*alpha)));
#if PRINT_UPDATES
  double tmp = (1 - lam)*u0 + lam*u1 + sh*sqrt(2 + lam*lam);
  printf("tri23(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, tmp);
  return tmp;
#else
  return (1 - lam)*u0 + lam*u1 + sh*sqrt(2 + lam*lam);
#endif
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tetra111(
  double u0, double u1, double u2,
  double s, double s0, double s1, double s2, double h) const
{
  (void) s0;
  (void) s1;
  (void) s2;

#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, h, s);
#endif

  double sh = s*h, du1 = u1 - u0, du2 = u2 - u0;
  double alpha1 = fabs(du1/sh), alpha2 = fabs(du2/sh);
  double A = alpha1*alpha1, B = alpha2*alpha2;
  double AQ1[6] = {2*A - 4, 2*A - 4, 2*A - 1, 4 - 2*A, 2 - 2*A, A - 1};
  double AQ2[6] = {2*B - 1, 2*B - 4, 2*B - 4, 2 - 2*B, 4 - 2*B, B - 1};

  int n = 0;
  double isects[8];
  for (int i = 0; i < 8; ++i) {
    isects[i] = -1;
  }
  this->intersect_conics(AQ1, AQ2, isects, n);

  double T = std::numeric_limits<double>::infinity();
  double lam0, lam1, lam2, l, lhs1, lhs2;
  for (int i = 0; i < 2*n; i += 2) {
    lam1 = isects[i];
    lam2 = isects[i + 1];
    if (lam1 < 0 || lam2 < 0 || lam1 + lam2 > 1) {
      continue;
    }
    lam0 = 1 - lam1 - lam2;
    l = sqrt(lam0*lam0 + lam1*lam1 + lam2*lam2);
    lhs1 = 2*lam1 + lam2 - 1 + du1*l/sh;
    lhs2 = lam1 + 2*lam2 - 1 + du2*l/sh;
    if (fabs(lhs1) < 1e-13 && fabs(lhs2) < 1e-13) {
      T = std::min(T, lam0*u0 + lam1*u1 + lam2*u2 + sh*l);
    }
  }
#if PRINT_UPDATES
  printf("tetra111(u0 = %g, u1 = %g, u2 = %g, s = %g, h = %g) -> %g\n",
         u0, u1, u2, s, h, T);
#endif
  return T;
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tetra122(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  using std::min;
  using std::sqrt;

  (void) s0;
  (void) s1;
  (void) s2;

#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, s, h);
#endif

  double T = std::numeric_limits<double>::infinity();

  double sh = s*h, du1 = u1 - u0, du2 = u2 - u0;
  double alpha1 = fabs(du1/sh), alpha1_sq = alpha1*alpha1;
  double alpha2 = fabs(du2/sh), alpha2_sq = alpha2*alpha2;

  if (alpha1_sq + alpha2_sq >= 1) {
    return T;
  }

  double denom = sqrt(1 - alpha1_sq - alpha2_sq);
  double lam1 = alpha1/denom;
  double lam2 = alpha2/denom;

  assert(lam1 >= 0);
  assert(lam2 >= 0);

  double lam0 = 1 - lam1 - lam2;
  if (lam0 < 0) {
    return T;
  }

  double l = sqrt(1 + lam1*lam1 + lam2*lam2);
  if (fabs(du1*l + sh*lam1) < 1e-13 && fabs(du2*l + sh*lam2) < 1e-13) {
    T = min(T, (1 - lam1 - lam2)*u0 + lam1*u1 + lam2*u2 + sh*l);
  }

#if PRINT_UPDATES
  printf("tetra122(u0 = %g, u1 = %g, u2 = %g, s = %g, h = %g) -> %g\n",
         u0, u1, u2, s, h, T);
#endif

  return T;
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tetra123(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  (void) s0;
  (void) s1;
  (void) s2;

#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, s, h);
#endif

  double sh = s*h, du1 = u1 - u0, du2 = u2 - u0;
  double beta1 = (du2 - 2*du1)/sh, beta1_sq = beta1*beta1;
  double beta2 = (du1 - du2)/sh, beta2_sq = beta2*beta2;
  double AQ1[6] = {beta1_sq - 1, 2*beta1_sq, 2*beta1_sq, 0, 0, beta1_sq};
  double AQ2[6] = {beta2_sq, 2*beta2_sq, 2*beta2_sq - 1, 0, 0, beta2_sq};

  int n = 0;
  double isects[8];
  this->intersect_conics(AQ1, AQ2, isects, n);

  double T = std::numeric_limits<double>::infinity();
  double lam1, lam2, l;
  for (int i = 0; i < n; i += 2) {
    lam1 = isects[i];
    lam2 = isects[i + 1];
    l = sqrt(lam1*lam1 + 2*lam2*(lam1 + lam2) + 1);
    if (fabs(lam1 + l*beta1) < 1e-13 && fabs(lam2 + l*beta2) < 1e-13) {
      T = std::min(T, (1 - lam1 - lam2)*u0 + lam1*u1 + lam2*u2 + sh*l);
    }
  }
#if PRINT_UPDATES
  printf("tetra123(u0 = %g, u1 = %g, u2 = %g, s = %g, h = %g) -> %g\n",
         u0, u1, u2, s, h, T);
#endif
  return T;
}

template <class rootfinder>
double olim3d_rhr_update_rules<rootfinder>::tetra222(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  (void) s0;
  (void) s1;
  (void) s2;

#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, s, h);
#endif

  double sh = s*h, du1 = u1 - u0, du2 = u2 - u0;
  double alpha1 = fabs(du1/sh), alpha2 = fabs(du2/sh);
  double A = alpha1*alpha1, B = alpha2*alpha2;
  double AQ1[6] = {2*A - 4, 2*A - 4, 2*A - 1, 4 - 2*A, 2 - 2*A, 2*A - 1};
  double AQ2[6] = {2*B - 1, 2*B - 4, 2*B - 4, 2 - 2*B, 4 - 2*B, 2*B - 1};

  int n = 0;
  double isects[8];
  for (int i = 0; i < 8; ++i) {
    isects[i] = -1;
  }
  this->intersect_conics(AQ1, AQ2, isects, n);

  double T = std::numeric_limits<double>::infinity();
  double lam1, lam2, l, lhs1, lhs2;
  for (int i = 0; i < 2*n; i += 2) {
    lam1 = isects[i];
    lam2 = isects[i + 1];
    if (lam1 < 0 || lam2 < 0 || lam1 + lam2 > 1) {
      continue;
    }
    l = sqrt(2*(lam1*lam1 + lam1*lam2 + lam2*lam2 - lam1 - lam2 + 1));
    lhs1 = 2*lam1 + lam2 - 1 + du1*l/sh;
    lhs2 = lam1 + 2*lam2 - 1 + du2*l/sh;
    if (fabs(lhs1) < 1e-13 && fabs(lhs2) < 1e-13) {
      T = std::min(T, (1 - lam1 - lam2)*u0 + lam1*u1 + lam2*u2 + sh*l);
    }
  }
#if PRINT_UPDATES
  printf("tetra222(u0 = %g, u1 = %g, u2 = %g, s = %g, h = %g) -> %g\n",
         u0, u1, u2, s, h, T);
#endif
  return T;
}

#endif // __OLIM_UPDATE_RULES_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
