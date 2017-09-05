#ifndef __OLIM26_RHR_IMPL_HPP__
#define __OLIM26_RHR_IMPL_HPP__

#include <algorithm>

template <class rootfinder>
double olim26_rhr_update_rules<rootfinder>::line1(
  double u0, double s, double s0, double h) const
{
  (void) s0;
  return u0 + s*h;
}

template <class rootfinder>
double olim26_rhr_update_rules<rootfinder>::line2(
  double u0, double s, double s0, double h) const
{
  (void) s0;
  return u0 + s*h*sqrt(2);
}

template <class rootfinder>
double olim26_rhr_update_rules<rootfinder>::line3(
  double u0, double s, double s0, double h) const
{
  (void) s0;
  return u0 + s*h*sqrt(3);
}

template <class rootfinder>
double olim26_rhr_update_rules<rootfinder>::tri12(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
  return rhr_diag(u0, u1, s, h);
}

template <class rootfinder>
double olim26_rhr_update_rules<rootfinder>::tri13(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
  double sh = s*h, alpha = fabs(u1 - u0)/sh, sgn = u0 > u1 ? 1 : -1;
  assert(2 > alpha*alpha);
  double lam = std::max(0.0, std::min(1.0, sgn*alpha/sqrt(2*(2 - alpha*alpha))));
  return (1 - lam)*u0 + lam*u1 + sh*sqrt(1 + 2*lam*lam);
}

template <class rootfinder>
double olim26_rhr_update_rules<rootfinder>::tri23(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
  double sh = s*h, alpha = fabs(u1 - u0)/sh, sgn = u0 > u1 ? 1 : -1;
  assert(1 > alpha*alpha);
  double lam = std::max(0.0, std::min(1.0, 2*sgn*alpha/sqrt(1 - alpha*alpha)));
  return (1 - lam)*u0 + lam*u1 + sh*sqrt(2 + lam*lam);
}

template <class rootfinder>
double olim26_rhr_update_rules<rootfinder>::tetra123(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  (void) s0;
  (void) s1;
  (void) s2;

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
  return T;
}

#endif

// __OLIM26_RHR_IMPL_HPP__
// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
