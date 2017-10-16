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
#include "common.macros.hpp"
#include "olim_rect_update_rules.hpp"
#include "olim_util.hpp"

double olim_rect_update_rules::line1(double u0, double s, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + s*h;
  printf("line1(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h;
#endif
}

double olim_rect_update_rules::line2(double u0, double s, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + s*h*sqrt2;
  printf("line2(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h*sqrt2;
#endif
}

double olim_rect_update_rules::line3(double u0, double s, double h)
  const
{
#if PRINT_UPDATES
  double tmp = u0 + s*h*sqrt3;
  printf("line3(u0 = %g, s = %g, h = %g) -> %g\n", u0, s, h, tmp);
  return tmp;
#else
  return u0 + s*h*sqrt3;
#endif
}

double olim_rect_update_rules::tri11(
  double u0, double u1, double s, double h) const
{
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, du = u1 - u0, alpha = fabs(du)/sh, alpha_sq = alpha*alpha;
  double T = INF(double);
  if (alpha_sq > 2) {
    return T;
  }
  double sgn = du < 0 ? 1 : -1;
  double lam = 0.5 + sgn*alpha/(2*sqrt(2 - alpha_sq));
  double l = sqrt(2*lam*(lam - 1) + 1);
  if (0 <= lam && lam <= 1) {
    T = (1 - lam)*u0 + lam*u1 + sh*l;
  }
#if PRINT_UPDATES
  printf("tri11(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
  return T;
}

double olim_rect_update_rules::tri12(double u0, double u1, double s, double h)
  const
{
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, du = u1 - u0, alpha = fabs(du)/sh, alpha_sq = alpha*alpha;
  double T = INF(double);
  if (alpha_sq > 1) {
    return T;
  }
  double lam = alpha/sqrt(1 - alpha_sq), l = sqrt(lam*lam + 1);
  if (0 <= lam && lam <= 1 && fabs(du*l + sh*lam) < 1e-13) {
    T = (1 - lam)*u0 + lam*u1 + sh*l;
  } else if (0 <= -lam && -lam <= 1 && fabs(du*l - sh*lam) < 1e-13) {
    T = (1 + lam)*u0 - lam*u1 + sh*l;
  }
#if PRINT_UPDATES
  printf("tri12(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
  return T;
}

double olim_rect_update_rules::tri13(double u0, double u1, double s, double h)
  const
{
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, du = u1 - u0, alpha = fabs(du)/sh, alpha_sq = alpha*alpha;
  double T = INF(double);
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

double olim_rect_update_rules::tri22(double u0, double u1, double s, double h)
  const
{
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, du = u1 - u0, alpha = fabs(du)/sh, alpha_sq = alpha*alpha;
  double T = INF(double);
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

double olim_rect_update_rules::tri23(double u0, double u1, double s, double h)
  const
{
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, h, s);
#endif
  double sh = s*h, du = u1 - u0, alpha = fabs(du)/sh, alpha_sq = alpha*alpha;
  double T = INF(double);
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

template <int M11, int M12, int M22, int E1, int E2, int F>
static double eval_Q(double lam1, double lam2) {
  return M11*lam1*lam1 + 2*M12*lam1*lam2 + M22*lam2*lam2 +
    2*E1*lam1 + 2*E2*lam2 + F;
}

/**
 * Compute the exact minimizing step size for the relaxed Newton's
 * method used to solve the minimization problem needed to do the
 * tetrahedral updates.
 */
template <int M11, int M12, int M22, int E1, int E2, int F>
static
double compute_step_size(
  double p1, double p2, double du1, double du2,
  double lam1, double lam2, double sh, bool & failure)
{
  failure = false;
  if (fabs(p1) <= EPS(double) && fabs(p2) <= EPS(double)) {
    return 0;
  }
  double du_dot_p = du1*p1 + du2*p2;
  double delta = std::pow(2*du_dot_p/sh, 2);
  double pMp = 2*(p1*(M11*p1 + M12*p2) + p2*(M12*p1 + M22*p2));
  double lMp = 2*(p1*(M11*lam1 + M12*lam2) + p2*(M12*lam1 + M22*lam2));
  double edotp = 2*(E1*p1 + E2*p2);
  double tmp1 = 2*pMp - delta, tmp2 = lMp + edotp;
  double a = pMp*tmp1/2;
  double b = tmp1*tmp2;
  double c = tmp2*tmp2 - delta*eval_Q<M11, M12, M22, E1, E2, F>(lam1, lam2);
  double disc = b*b - 4*a*c;
  if (disc < 0) {
    failure = true;
    return INF(double);
  }
  double lhs = -b/(2*a);
  double rhs = sqrt(disc)/(2*a);
  return lhs + rhs <= 0 ? lhs + rhs : lhs - rhs;
}

using compute_p_func =
  void(*)(double, double, double, double, double, double &, double &);

template <int M11, int M12, int M22, int E1, int E2, int F>
static
double tetra_newton(
  double lam1init, double lam2init, double u0, double u1, double u2, double sh,
  compute_p_func compute_p)
{
  using std::max;
  double du1 = u1 - u0, du2 = u2 - u0;
  double lam1 = lam1init, lam2 = lam2init;
  double p1, p2, Q, relerr, alpha;

  /**
   * Do a relaxed Newton iteration to minimize the line integral
   * required by the tetrahedral update.
   */
  int niters = 0;
  bool failure = false;
  do {
    /**
     * Compute the iteration step, including exact step size.
     */
    Q = eval_Q<M11, M12, M22, E1, E2, F>(lam1, lam2);
    compute_p(lam1, lam2, du1, du2, sqrt(Q)/sh, p1, p2);
    alpha = compute_step_size<M11, M12, M22, E1, E2, F>(
      p1, p2, du1, du2, lam1, lam2, sh, failure);
    if (failure) return INF(double);
    if (fabs(alpha) <= EPS(double)) break;

    /**
     * Update lambda by the computed step size and abort if lambda
     * lies outside of the 2-simplex.
     */
    lam1 += alpha*p1 , lam2 += alpha*p2;
    if (lam1 < 0 || lam2 < 0 || 1 - lam1 - lam2 < 0) return INF(double);

    /**
     * Compute the error for the iteration.
     */
    relerr = max(fabs(p1), fabs(p2))/max(fabs(lam1), fabs(lam2));

    /**
     * Keep track of the iteration number and print some debug
     * information if necessary.
     */
    ++niters;
    if (niters > 10) {
#ifdef EIKONAL_DEBUG
      printf("u0 = %0.16g, u1 = %0.16g, u2 = %0.16g, sh = %0.16g\n",
             u0, u1, u2, sh);
      std::abort();
#else
      break;
#endif
    }
  } while (relerr > 1e-15);

  /**
   * Compute the minimized quantity and return it.
   */
  Q = eval_Q<M11, M12, M22, E1, E2, F>(lam1, lam2);
  return (1 - lam1 - lam2)*u0 + lam1*u1 + lam2*u2 + sh*sqrt(Q);
}

static void p111(double lam1, double lam2, double du1, double du2,
                 double delta, double & p1, double & p2) {
  double H11 = delta*(1 + lam1*(3*lam1 - 2));
  double H12 = delta*(lam1*(3*lam2 - 1) - lam2);
  double H22 = delta*(1 + lam2*(3*lam2 - 2));
  double G1 = du1 + (2*lam1 + lam2 - 1)/delta;
  double G2 = du2 + (lam1 + 2*lam2 - 1)/delta;
  p1 = H11*G1 + H12*G2;
  p2 = H12*G1 + H22*G2;
}

double olim_rect_update_rules::tetra111(
  double u0, double u1, double u2, double s, double h) const
{
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, h, s);
#endif
  double T = tetra_newton<2, 1, 2, -1, -1, 1>(
    1.0/3.0, // lam1init
    1.0/3.0, // lam1init
    u0,      // u0
    u1,      // u1
    u2,      // u2
    s*h,     // sh
    p111);   // compute_p
#if PRINT_UPDATES
  printf("tetra111(u0 = %0.16g, u1 = %0.16g, u2 = %0.16g, s = %0.16g, "
         "h = %0.16g) -> %g\n", u0, u1, u2, s, h, T);
#endif
  return T;
}

double olim_rect_update_rules::tetra122(
  double u0, double u1, double u2, double s, double h) const
{
  using std::min;
  using std::sqrt;

#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, s, h);
#endif

  double T = INF(double);

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

static void p123(double lam1, double lam2, double du1, double du2,
                 double delta, double & p1, double & p2) {
  double H11 = delta*(2 + lam1*lam1);
  double H12 = delta*(lam1*lam2 - 1);
  double H22 = delta*(1 + lam2*lam2);
  double G1 = du1 + (lam1 + lam2)/delta;
  double G2 = du2 + (lam1 + 2*lam2)/delta;
  p1 = H11*G1 + H12*G2;
  p2 = H12*G1 + H22*G2;
}

double olim_rect_update_rules::tetra123(
  double u0, double u1, double u2, double s, double h) const
{
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, s, h);
#endif
  double T = tetra_newton<1, 1, 2, 0, 0, 1>(
    0.0,   // lam1init
    0.0,   // lam2init
    u0,    // u0
    u1,    // u1
    u2,    // u2
    s*h,   // sh
    p123); // compute_p
#if PRINT_UPDATES
  printf("tetra123(u0 = %g, u1 = %g, u2 = %g, s = %g, h = %g) -> %g\n",
         u0, u1, u2, s, h, T);
#endif
  return T;
}

static void p222(double lam1, double lam2, double du1, double du2,
                 double delta, double & p1, double & p2) {
  double H11 = delta*(3 + lam1*(3*lam1 - 2))/4.0;
  double H12 = delta*(lam1*(3*lam2 - 1) - lam2 - 1)/4.0;
  double H22 = delta*(3 + lam2*(3*lam2 - 2))/4.0;
  double G1 = du1 + (2*lam1 + lam2 - 1)/delta;
  double G2 = du2 + (lam1 + 2*lam2 - 1)/delta;
  p1 = H11*G1 + H12*G2;
  p2 = H12*G1 + H22*G2;
}

double olim_rect_update_rules::tetra222(
  double u0, double u1, double u2, double s, double h) const
{
#ifdef EIKONAL_DEBUG
  check_params(u0, u1, u2, s, h);
#endif
  double T = tetra_newton<2, 1, 2, -1, -1, 2>(
    1.0/3.0, // lam1init
    1.0/3.0, // lam2init
    u0,      // u0
    u1,      // u1
    u2,      // u2
    s*h,     // sh
    p222);   // compute_p
#if PRINT_UPDATES
  printf("tetra123(u0 = %g, u1 = %g, u2 = %g, s = %g, h = %g) -> %g\n",
         u0, u1, u2, s, h, T);
#endif
  return T;
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
