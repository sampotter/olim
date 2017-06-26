#include "olim_8pt_util.hpp"

#include <cassert>
#include <cmath>
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
