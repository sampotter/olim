#include "olim6_rhr.hpp"

#include "olim_util.hpp"

template <class rootfinder>
double olim6_rhr_update_rules<rootfinder>::adj1pt(
  double u0, double s, double s0, double h) const
{
  (void) s0;
  return u0 + s*h;
}

template <class rootfinder>
double olim6_rhr_update_rules<rootfinder>::adj2pt(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
  return rhr_adj(u0, u1, s, h);
}

template <class rootfinder>
double olim6_rhr_update_rules<rootfinder>::adj3pt(
  double u0, double u1, double u2,
  double s, double s0, double s1, double s2, double h) const
{
  (void) s0;
  (void) s1;
  (void) s2;

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
  return T;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
