#ifndef __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#if PRINT_UPDATES
#    include <cstdio>
#endif
#include <cmath>

#include "common.hpp"
#include "common.defs.hpp"
#include "common.macros.hpp"
#include "olim_util.hpp"

template <class speed_est, char degree>
template <char p0, char p1>
double
update_rules::tri_updates<speed_est, degree>::tri(
  double u0, double u1, double s, double s0, double s1, double h,
  ffvec<p0>, ffvec<p1>, double tol) const
{
  double T = tri_impl(
    u0, u1, s, s0, s1, h, ffvec<p0> {}, ffvec<p1> {}, tol,
    std::integral_constant<char, degree> {});
#if PRINT_UPDATES
  printf("tri<%d, %d, %d>::update_impl(u0 = %g, u1 = %g, "
         "s = %g, s0 = %g, s1 = %g, h = %g) -> %g\n", degree, p0, p1,
         u0, u1, s, s0, s1, h, T);
#endif
  return T;
}

#define l__(x) std::sqrt((dp_dot_dp*(x) + 2*dp_dot_p0)*(x) + p0_dot_p0)
#define check__(x) std::fabs(alpha*l__(x) - dp_dot_p0 - dp_dot_dp*(x))
#define char_sqrt__(c) _sqrt_table[static_cast<int>(c)]
#define l0__(i) (u##i + sh*char_sqrt__(p##i##_dot_p##i))

/**
 * F0 specialization
 */
template <class speed_est, char degree>
template <char p0, char p1>
double
update_rules::tri_updates<speed_est, degree>::tri_impl(
  double u0, double u1, double s, double s0, double s1, double h,
  ffvec<p0>, ffvec<p1>, double tol, std::integral_constant<char, 0>) const
{
  using std::min;

  (void) tol;

  static constexpr double _sqrt_table[4] = {
    0.0,
    1.0,
    1.4142135623730951,
    1.7320508075688772
  };

  constexpr char p0_dot_p0 = dot(p0, p0);
  constexpr char p0_dot_p1 = dot(p0, p1);
  constexpr char p1_dot_p1 = dot(p1, p1);
  constexpr char dp_dot_p0 = p0_dot_p1 - p0_dot_p0;
  constexpr char dp_dot_p0_sq = dp_dot_p0*dp_dot_p0;
  constexpr char dp_dot_dp = p1_dot_p1 - 2*p0_dot_p1 + p0_dot_p0;

  double const du = u1 - u0;
  double const sh = this->s_hat(s, s0, s1)*h;
  double const alpha = -du/sh, alpha_sq = alpha*alpha;
  double const tmp = alpha_sq - dp_dot_dp;
  double const a = dp_dot_dp*tmp;
  double const b = dp_dot_p0*tmp;
  double const c = alpha_sq*p0_dot_p0 - dp_dot_p0_sq;
  double const disc = b*b - a*c;

  if (disc < 0 || a == 0) {
    return std::min(l0__(0), l0__(1));
  } else {
    double const lhs = -b/a, rhs = sqrt(disc)/a;
    double const lam1 = lhs - rhs, lam2 = lhs + rhs;
    double const lam = check__(lam1) < check__(lam2) ? lam1 : lam2;
    return lam < 0 || 1 < lam ?
      min(l0__(0), l0__(1)) :
      u0 + lam*du + sh*l__(lam);
  }
}

#undef l__
#undef check__
#undef char_sqrt__
#undef l0__

#define u__(x) ((1 - (x))*u0 + (x)*u1)
#define q__(x) ((dp_dot_dp*x + 2*dp_dot_p0)*(x) + p0_dot_p0)
#define l__(x) std::sqrt(q__(x))
#define s__(x) ((1 - theta)*s + theta*((1 - (x))*s0 + (x)*s1))
#define F1__(x) (u__(x) + h*s__(x)*l__(x))
#define dF1__(x) (du + h*(ds*theta*q__(x) + s__(x)*dp_dot_plam)/l__(x))
#define d2F1__(x) h*(s__(x)*dp_dot_plam*dp_dot_plam/q__(x) + \
                     ds*theta*dp_dot_plam + 2*s__(x)*dp_dot_dp)/l__(x)

/**
 * F1 specialization
 */
template <class speed_est, char degree>
template <char p0, char p1>
double
update_rules::tri_updates<speed_est, degree>::tri_impl(
  double u0, double u1, double s, double s0, double s1, double h,
  ffvec<p0>, ffvec<p1>, double tol, std::integral_constant<char, 1>) const
{
  constexpr char p0_dot_p0 = dot(p0, p0);
  constexpr char p0_dot_p1 = dot(p0, p1);
  constexpr char p1_dot_p1 = dot(p1, p1);
  constexpr char dp_dot_p0 = p0_dot_p1 - p0_dot_p0;
  constexpr char dp_dot_dp = p1_dot_p1 - 2*p0_dot_p1 + p0_dot_p0;

  constexpr double c1 = 1e-4;

  double const ds = s1 - s0, du = u1 - u0, theta = this->theta();

  bool conv;
  double lam[2], F1[2], dF1, d2F1, g, dp_dot_plam, alpha;
  lam[0] = 0.5;
  F1[0] = F1__(lam[0]);
  do {
    alpha = 1;
    dp_dot_plam = dp_dot_p0 + lam[0]*dp_dot_dp;
    dF1 = dF1__(lam[0]);
    d2F1 = d2F1__(lam[0]);
    g = -dF1/d2F1;

    double lam_ = lam[0] + alpha*g;
    double u_ = u__(lam_);
    double s_ = s__(lam_);
    double l_ = l__(lam_);
    double tmp = u_ + h*s_*l_;
    while (tmp > F1[0] + c1*alpha*dF1*g) {
      alpha *= 0.9;
      lam_ = lam[0] + alpha*g;
      u_ = u__(lam_);
      s_ = s__(lam_);
      l_ = l__(lam_);
      tmp = u_ + h*s_*l_;
    }
    lam[1] = std::max(0., std::min(1., lam[0] + alpha*g));
    F1[1] = F1__(lam[1]);
    conv = fabs(lam[1] - lam[0]) <= tol || fabs(F1[1] - F1[0]) <= tol;
    lam[0] = lam[1];
    F1[0] = F1[1];
  } while (!conv);
  return F1[1];
}

#undef u__
#undef q__
#undef l__
#undef theta__
#undef s__
#undef F1__
#undef dF1__
#undef d2F1__

#endif // __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__
