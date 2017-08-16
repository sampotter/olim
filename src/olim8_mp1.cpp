#include "olim8_mp1.hpp"

#include <cmath>
#include <cstdio>

#include "olim_util.hpp"
#include "rootfinding.hpp"

double olim8_mp1_update_rules::adj1pt(double u0, double s, double s0,
                                      double h) const {
  return u0 + h*(s + s0)/2;
}

double olim8_mp1_update_rules::diag1pt(double u0, double s, double s0,
                                       double h) const {
  return u0 + h*(s + s0)*std::sqrt(2)/2;
}

double olim8_mp1_gsl_update_rules::adj2pt(double u0, double u1, double s,
                                          double s0, double s1,
                                          double h) const {
  double sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
  return sbar0 == sbar1 ? rhr_adj(u0, u1, sbar0, h) :
    mp1_adj(u0, u1, sbar0, sbar1, h);
}

double olim8_mp1_gsl_update_rules::diag2pt(double u0, double u1, double s,
                                           double s0, double s1,
                                           double h) const {
  double sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
  return sbar0 == sbar1 ? rhr_diag(u0, u1, sbar0, h) :
    mp1_diag(u0, u1, sbar0, sbar1, h);
}

static double evalquartic(double x, double * a) {
  return a[0] + x*(a[1] + x*(a[2] + x*(a[3] + a[4]*x)));
}

double olim8_mp1_bsearch_update_rules::adj2pt(double u0, double u1, double s,
                                              double s0, double s1,
                                              double h) const {
  double sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
  if (sbar0 == sbar1) {
    return rhr_adj(u0, u1, sbar0, h);
  }

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

  double lam, Tnew, argmin, roots[5], T = std::numeric_limits<double>::infinity();
  find_quartic_roots(a, roots);
  (void) argmin;

  int i = 0;
  while ((lam = roots[i++]) != -1) {
    double lhs = (u0 - u1)*std::sqrt(1 - 2*lam + 2*lam*lam)/h;
    double rhs = -s0*(4*lam*lam - 5*lam + 2) + s1*(4*lam*lam - 3*lam + 1);
    if (fabs(lhs - rhs)/fabs(lhs) < 1e-6) {
      Tnew = (1 - lam)*u0+ lam*u1 +
        ((1 - lam)*s0 + lam*s1)*h*std::sqrt(1 - 2*lam + 2*lam*lam + 1);
      if (Tnew < T) {
        T = Tnew;
        argmin = lam;
      }
    }
  }
  // if (std::isinf(T)) {
  //   printf("adj2pt: T not updated\n");
  // } else {
  //   printf("adj2pt: T <- %g (argmin = %g)\n", T, argmin);
  // }
  return T;
}

double olim8_mp1_bsearch_update_rules::diag2pt(double u0, double u1, double s,
                                               double s0, double s1,
                                               double h) const {
  double sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
  if (sbar0 == sbar1) {
    return rhr_diag(u0, u1, sbar0, h);
  }

  double alpha = std::fabs((u0 - u1)/h), alpha_sq = alpha*alpha;
  double ds = s1 - s0, ds_sq = ds*ds;
  double a[] = {
    ds_sq - alpha_sq,
    2*s0*ds,
    4*ds_sq + s0*s0 - alpha_sq,
    4*s0*ds,
    4*ds_sq
  };

  double lam, Tnew, argmin, roots[5],
    T = std::numeric_limits<double>::infinity();
  find_quartic_roots(a, roots);
  (void) argmin;

  int i = 0;
  while ((lam = roots[i++]) != -1) {
    double lhs = (u0 - u1)*std::sqrt(lam*lam + 1)/h;
    double rhs = s1 + 2*s1*lam*lam - s0*(1 - lam + 2*lam*lam);
    if (fabs(lhs - rhs)/fabs(lhs) < 1e-6) {
      Tnew = (1 - lam)*u0+ lam*u1 +
        ((1 - lam)*s0 + lam*s1)*h*std::sqrt(lam*lam + 1);
      if (Tnew < T) {
        T = Tnew;
        argmin = lam;
      }
    }
  }
  // if (std::isinf(T)) {
  //   printf("diag2pt: T not updated\n");
  // } else {
  //   printf("diag2pt: T <- %g (argmin = %g)\n", T, argmin);
  // }
  return T;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
