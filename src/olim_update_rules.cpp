#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.defs.hpp"
#include "olim_update_rules.hpp"
#include "olim_util.hpp"

double olim3d_rhr_update_rules::line1(double u0, double s, double s0, double h)
  const
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

double olim3d_rhr_update_rules::line2(double u0, double s, double s0, double h)
  const
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

double olim3d_rhr_update_rules::line3(double u0, double s, double s0, double h)
  const
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

// TODO: make this solve the unconstrained problem
double olim3d_rhr_update_rules::tri11(
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

// TODO: make this solve the unconstrained problem
double olim3d_rhr_update_rules::tri12(
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

double olim3d_rhr_update_rules::tri13(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, du = u1 - u0, alpha = fabs(du)/sh, alpha_sq = alpha*alpha;
  double T = std::numeric_limits<double>::infinity();
  if (alpha_sq > 2) {
    return T;
  }
  double lam = alpha/sqrt(2*(2 - alpha_sq)), l = sqrt(1 + 2*lam*lam);
  if (0 <= lam && lam <= 1 && fabs(du*l + sh*lam) < 1e-13) {
    T = (1 - lam)*u0 + lam*u1 + sh*l;
  } else if (0 <= -lam && -lam <= 1 && fabs(du*l - sh*lam) < 1e-13) {
    T = (1 + lam)*u0 - lam*u1 + sh*l;
  }
#if PRINT_UPDATES
  printf("tri13(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
  return T;
}

double olim3d_rhr_update_rules::tri22(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, du = u1 - u0, alpha = fabs(du)/sh, alpha_sq = alpha*alpha;
  double T = std::numeric_limits<double>::infinity();
  if (alpha_sq > 2) {
    return T;
  }
  double rhs = sqrt3*alpha/(2*sqrt(4 - alpha_sq));
  double lam = 0.5 - rhs, l = sqrt(2*(1 - lam*(1 - lam)));
  if (0 <= lam && lam <= 1 && fabs(du*l + sh*lam) < 1e-13) {
    T = (1 - lam)*u0 + lam*u1 + sh*l;
  } else {
    lam = 0.5 + rhs, l = sqrt(2*(1 - lam*(1 - lam)));
    if (0 <= lam && lam <= 1 && fabs(du*l * sh*lam) < 1e-13) {
      T = (1 - lam)*u0 + lam*u1 + sh*l;
    }
  }
#if PRINT_UPDATES
  printf("tri22(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
  return T;
}

double olim3d_rhr_update_rules::tri23(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, du = u1 - u0, alpha = fabs(du)/sh, alpha_sq = alpha*alpha;
  double T = std::numeric_limits<double>::infinity();
  if (alpha_sq > 1) {
    return T;
  }
  double lam = sqrt2*alpha/sqrt(1 - alpha_sq), l = sqrt(2 + lam*lam);
  if (0 <= lam && lam <= 1 && fabs(du*l + sh*lam) < 1e-13) {
    T = (1 - lam)*u0 + lam*u1 + sh*l;
  } else if (0 <= -lam && -lam <= 1 && fabs(du*l - sh*lam) < 1e-13) {
    T = (1 + lam)*u0 - lam*u1 + sh*l;
  }
#if PRINT_UPDATES
  printf("tri23(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
  return T;
}

double olim3d_rhr_update_rules::tetra111(
  double u0, double u1, double u2,
  double s, double s0, double s1, double s2, double h) const
{
  using std::max;

  (void) s0;
  (void) s1;
  (void) s2;

#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, h, s);
#endif

  double sh = s*h, du1 = u1 - u0, du2 = u2 - u0;
  double T = std::numeric_limits<double>::infinity();
  double x = 1./3., y = 1./3., p1, p2, err, alpha;

  int niters = 0;
  do {
    double l = sqrt(pow(1 - x - y, 2) + x*x + y*y);
    double c = l/sh;
    double H11 = c*(1 + x*(3*x - 2)), H12 = c*(x*(3*y - 1) - y),
      H22 = c*(1 + y*(3*y - 2));
    double G1 = du1 + (2*x + y - 1)/c, G2 = du2 + (x + 2*y - 1)/c;

    p1 = H11*G1 + H12*G2;
    p2 = H12*G1 + H22*G2;

    // compute step size---clean this up later
    if (fabs(p1) <= std::numeric_limits<double>::epsilon() &&
        fabs(p2) <= std::numeric_limits<double>::epsilon()) {
      break;
    } else {
      double du_dot_p = du1*p1 + du2*p2;
      double delta = 2*du_dot_p/sh;
      delta *= delta;

      double pMp = p1*(4*p1 + 2*p2) + p2*(2*p1 + 4*p2);
      double lMp = p1*(4*x + 2*y) + p2*(2*x + 4*y);
      double edotp = -2*(p1 + p2);

      double tmp1 = 2*pMp - delta, tmp2 = lMp + edotp, tmp3 = 1 - x - y;
      double a = pMp*tmp1/2;
      double b = tmp1*tmp2;
      double c = tmp2*tmp2 - delta*(tmp3*tmp3 + x*x + y*y);
      double disc = b*b - 4*a*c;

      assert(disc >= 0);
      double lhs = -b/(2*a);
      double rhs = sqrt(disc)/(2*a);
      alpha = lhs + rhs <= 0 ? lhs + rhs : lhs - rhs;
    }
    if (fabs(alpha) <= std::numeric_limits<double>::epsilon()) {
      break;
    }

    x += alpha*p1 , y += alpha*p2;
    if (x < 0 || y < 0 || 1 - x - y < 0) {
      goto coda;
    }
    err = max(fabs(p1), fabs(p2))/max(fabs(x), fabs(y));
    ++niters;
    if (niters > 10) {
      break;
      // printf("u0 = %0.16g, u1 = %0.16g, u2 = %0.16g, s = %0.16g, h = %0.16g\n",
      //        u0, u1, u2, s, h);
      // std::abort();
    }
  } while (err > 1e-15);
  T = (1 - x - y)*u0 + x*u1 + y*u2 + sh*sqrt(pow(1 - x - y, 2) + x*x + y*y);

  coda:
#if PRINT_UPDATES
  printf("tetra111(u0 = %0.16g, u1 = %0.16g, u2 = %0.16g, s = %0.16g, "
         "h = %0.16g) -> %g\n", u0, u1, u2, s, h, T);
#endif
  return T;
}

double olim3d_rhr_update_rules::tetra122(
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

double olim3d_rhr_update_rules::tetra123(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  using std::max;

  (void) s0;
  (void) s1;
  (void) s2;

#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, s, h);
#endif

  double sh = s*h, du1 = u1 - u0, du2 = u2 - u0;

  double T = std::numeric_limits<double>::infinity();
  double x = 0, y = 0, dx, dy, err;
  do {
    double l = sqrt(x*x + 2*x*y + 2*y*y + 1);
    double c = l/sh;
    double H11 = c*(2 + x*x), H12 = c*(x*y - 1), H22 = c*(1 + y*y);
    double G1 = du1 + (x + y)/c, G2 = du2 + (x + 2*y)/c;

    dx = -(H11*G1 + H12*G2); x += dx;
    dy = -(H12*G1 + H22*G2), y += dy;
    if (x < 0 || y < 0 || 1 - x - y < 0) {
      goto coda;
    }
    err = max(fabs(dx), fabs(dy))/max(fabs(x), fabs(y));
  } while (err > 1e-15);
  T = (1 - x - y)*u0 + x*u1 + y*u2 + sh*sqrt(x*x + 2*x*y + 2*y*y + 1);

  coda:
#if PRINT_UPDATES
  printf("tetra123(u0 = %g, u1 = %g, u2 = %g, s = %g, h = %g) -> %g\n",
         u0, u1, u2, s, h, T);
#endif
  return T;
}

double olim3d_rhr_update_rules::tetra222(
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  using std::max;

  (void) s0;
  (void) s1;
  (void) s2;

#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, s, h);
#endif

  double sh = s*h, du1 = u1 - u0, du2 = u2 - u0;
  double T = std::numeric_limits<double>::infinity();
  double x = 1./3., y = 1./3., dx, dy, err;
  do {
    double l = sqrt(2*(1 + x*x + x*(y - 1) - y + y*y));
    double c = l/sh;
    double H11 = c*(3 + x*(3*x - 2))/4.0, H12 = c*(x*(3*y - 1) - y - 1)/4.0,
      H22 = c*(3 + y*(3*y - 2))/4.0;
    double G1 = du1 + (2*x + y - 1)/c, G2 = du2 + (x + 2*y - 1)/c;

    dx = -(H11*G1 + H12*G2); x += dx;
    dy = -(H12*G1 + H22*G2), y += dy;
    if (x < 0 || y < 0 || 1 - x - y < 0) {
      goto coda;
    }
    err = max(fabs(dx), fabs(dy))/max(fabs(x), fabs(y));
  } while (err > 1e-15);
  T = (1 - x - y)*u0 + x*u1 + y*u2 + sh*sqrt(2*(1 + x*x + x*(y - 1) - y + y*y));

  coda:
#if PRINT_UPDATES
  printf("tetra222(u0 = %g, u1 = %g, u2 = %g, s = %g, h = %g) -> %g\n",
         u0, u1, u2, s, h, T);
#endif
  return T;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
