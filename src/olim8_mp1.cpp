#include "olim8_mp1.hpp"

#include "common.defs.hpp"
#include "olim_util.hpp"
#include "qroots.hpp"

double olim8_mp1_update_rules::tri11(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  double sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
  if (sbar0 == sbar1) {
    return rhr_adj(u0, u1, sbar0, h);
  }

  double alpha_sq = std::pow((u0 - u1)/h, 2);
  double sbar0_sq = sbar0*sbar0;
  double dsbar = sbar1 - sbar0;
  double dsbar_sq = dsbar*dsbar;
  double a[] = {
    sbar0_sq - alpha_sq - 2*sbar0*dsbar + dsbar_sq,
    -4*sbar0_sq + 2*alpha_sq + 10*sbar0*dsbar - 6*dsbar_sq,
    4*sbar0_sq - 2*alpha_sq - 20*sbar0*dsbar + 17*dsbar_sq,
    16*sbar0*dsbar - 24*dsbar_sq,
    16*dsbar_sq
  };

  double lam, Tnew, argmin, roots[4] = {-1, -1, -1, -1},
    T = std::numeric_limits<double>::infinity(), lhs, rhs;
  qroots(a, roots);

  int i = 0;
  while ((lam = roots[i++]) != -1) {
    lhs = (u0 - u1)*std::sqrt(1 - 2*lam + 2*lam*lam)/h;
    rhs = -sbar0*(4*lam*lam - 5*lam + 2) + sbar1*(4*lam*lam - 3*lam + 1);
    if (fabs(lhs - rhs)/fabs(lhs) < 1e-6) {
      Tnew = (1 - lam)*u0+ lam*u1 +
        ((1 - lam)*sbar0 + lam*sbar1)*h*std::sqrt(1 - 2*lam + 2*lam*lam + 1);
      if (Tnew < T) {
        T = Tnew;
        argmin = lam;
      }
    }
  }

  (void) argmin;

  return T;
}

double olim8_mp1_update_rules::tri12(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  double sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
  if (sbar0 == sbar1) {
    return rhr_diag(u0, u1, sbar0, h);
  }

  double alpha = std::fabs((u0 - u1)/h);
  double dsbar = sbar1 - sbar0;
  double const a[] = {
    (dsbar - alpha)*(dsbar + alpha),
    2*sbar0*dsbar,
    4*dsbar*dsbar + (sbar0 - alpha)*(sbar0 + alpha),
    4*sbar0*dsbar,
    4*dsbar*dsbar
  };

  double lam, Tnew, argmin, roots[4] = {-1, -1, -1, -1},
    T = std::numeric_limits<double>::infinity(), lhs, rhs;
  qroots(a, roots);

  int i = 0;
  while ((lam = roots[i++]) != -1) {
    lhs = (u0 - u1)*std::sqrt(lam*lam + 1)/h;
    rhs = sbar1 + 2*sbar1*lam*lam - sbar0*(1 - lam + 2*lam*lam);
    if (fabs(lhs - rhs)/fabs(lhs) < 1e-6) {
      Tnew = (1 - lam)*u0+ lam*u1 +
        ((1 - lam)*sbar0 + lam*sbar1)*h*std::sqrt(lam*lam + 1);
      if (Tnew < T) {
        T = Tnew;
        argmin = lam;
      }
    }
  }
  
  (void) argmin;
  
  return T;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
