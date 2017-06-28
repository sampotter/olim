#include "olim8_util.hpp"

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

double rhr_adj(double u0, double u1, double s_est, double h) {
  check_params(u0, u1, h);
  
  double c = s_est*std::fabs(u0 - u1)/h;
  if (c > std::sqrt(2)/2) {
    return std::numeric_limits<double>::infinity();
  }
  double rad = sqrt(3)*c/sqrt(c*c + 2);
  double lam1 = (1 + rad)/2;
  double lam2 = (1 - rad)/2;
  assert(std::fabs(1 - lam1 - lam2) < 1e-15 ||
         (0 < lam1 && lam1 < 1) != (0 < lam2 && lam2 < 1));
  double lam = 0 < lam1 && lam1 < 1 ? lam1 : lam2;
  return (1 - lam)*u0 + lam*u1 + s_est*h*sqrt(2*lam*(1 - lam) + 1);
}

double rhr_diag(double u0, double u1, double s_est, double h) {
  check_params(u0, u1, h);
  assert(s_est >= 0);
  
  double c = s_est*std::fabs(u0 - u1)/h;
  if (c > std::sqrt(2)/2) {
    return std::numeric_limits<double>::infinity();
  }
  double lam = c/std::sqrt(1 - c*c);
  assert(0 <= lam && lam <= 1);
  return (1 - lam)*u0 + lam*u1 + s_est*h*sqrt(lam*lam + 1);
}

double mp1_adj(double u0, double u1, double sbar0, double sbar1, double h) {
  check_params(u0, u1, h);
  assert(sbar0 >= 0);
  assert(sbar1 >= 0);

  double alpha_sq = std::pow((u0 - u1)/h, 2);
  double dsbar = sbar1 - sbar0;
  double dsbar_sq = dsbar*dsbar;
  double sbar0_sq = sbar0*sbar0;
  double sbar0_times_dsbar = sbar0*dsbar;
  double a[] = {
    sbar0_sq - alpha_sq - 2*sbar0_times_dsbar,
    -4*sbar0_sq + 2*alpha_sq + 10*sbar0_times_dsbar - 6*dsbar_sq,
    4*sbar0_sq - 2*alpha_sq - 20*sbar0_times_dsbar + 17*dsbar_sq,
    16*sbar0_times_dsbar - 24*dsbar_sq,
    16*dsbar_sq
  };
  // printf("mp1_adj: {");
  // for (int i = 0; i < 4; ++i) printf("%g, ", a[i]);
  // printf("%g}\n", a[4]);

  double z[8];
  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(5);
  gsl_poly_complex_solve(a, 5, w, z);
  gsl_poly_complex_workspace_free(w);
  // printf("ROOTS (u0 = %g, u1 = %g, sbar0 = %g, sbar1 = %g):\n", u0, u1, sbar0, sbar1);
  // for (int i = 0; i < 5; i++) {
  //   printf("z%d = %+.18f %+.18f\n", i, z[2*i], z[2*i + 1]);
  // }

  double lam = -1;
  for (int i = 0; i < 4; ++i) {
    if (z[2*i + 1] != 0) continue;
    lam = z[2*i];
    if (0 <= lam && lam <= 1) break;
  }
  if (0 <= lam && lam <= 1) {
    return (1 - lam)*u0 + lam*u1 +
      h*((1 - lam)*sbar0 + lam*sbar1)*std::sqrt(2*lam*lam - 2*lam + 1);
  } else {
    return rhr_adj(u0, u1, (sbar0 + sbar1)/2, h);
  }
}

double mp1_diag(double u0, double u1, double sbar0, double sbar1, double h) {
  check_params(u0, u1, h);
  assert(sbar0 >= 0);
  assert(sbar1 >= 0);

  double alpha_sq = std::pow((u0 - u1)/h, 2);
  double dsbar = sbar1 - sbar0;
  double a[] = {
    dsbar - alpha_sq,
    2*sbar0*dsbar,
    4*dsbar + sbar0*sbar0 - alpha_sq,
    4*sbar0*dsbar,
    4*dsbar
  };
  // printf("mp1_diag: {");
  // for (int i = 0; i < 4; ++i) printf("%g, ", a[i]);
  // printf("%g}\n", a[4]);

  double z[8];
  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(5);
  gsl_poly_complex_solve(a, 5, w, z);
  gsl_poly_complex_workspace_free(w);
  // printf("ROOTS (u0 = %g, u1 = %g, sbar0 = %g, sbar1 = %g):\n", u0, u1, sbar0, sbar1);
  // for (int i = 0; i < 5; i++) {
  //   printf("z%d = %+.18f %+.18f\n", i, z[2*i], z[2*i + 1]);
  // }

  double lam = -1;
  for (int i = 0; i < 4; ++i) {
    if (z[2*i + 1] != 0) continue;
    lam = z[2*i];
    if (0 <= lam && lam <= 1) break;
  }
  if (0 <= lam && lam <= 1) {
    return (1 - lam)*u0 + lam*u1 +
      h*((1 - lam)*sbar0 + lam*sbar1)*std::sqrt(1 + lam*lam);
  } else {
    return rhr_diag(u0, u1, (sbar0 + sbar1)/2, h);
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
