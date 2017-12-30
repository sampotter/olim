#ifndef __UPDATE_RULES_TETRA_UPDATES_UTIL_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_UTIL_HPP__

#include "update_rules.utils.hpp"

template <char p0, char p1, char p2>
constexpr char dpi_dot_p0(char i) {
  assert(i == 1 || i == 2);
  if (i == 1) {
    return dot(p1, p0) - dot(p0, p0);
  } else {
    return dot(p2, p0) - dot(p0, p0);
  }
}

/**
 * Coefficients of conic section for d = 2.
 */
template <char p0, char p1, char p2>
constexpr char coef(char i) {
  constexpr char lut[6] = {
    dot(p1, p1) - 2*dot(p0, p1) + dot(p0, p0),
    dot(p1, p2) - dot(p1, p0) - dot(p2, p0) + dot(p0, p0),
    dot(p2, p2) - 2*dot(p0, p2) + dot(p0, p0),
    dpi_dot_p0<p0, p1, p2>(1),
    dpi_dot_p0<p0, p1, p2>(2),
    dot(p0, p0)
  };
  assert(0 <= i && i < 6);
  return lut[static_cast<int>(i)];
}

/**
 * Coefficient of the $2 \times 2$ matrix \left({\delta P}^\top
 * {\delta P}\right)^{-1}.
 */
template <char p0, char p1, char p2>
constexpr double dPt_times_dP_coef(char i) {
  assert(i == 0 || i == 1 || i == 2);
  constexpr char a = coef<p0, p1, p2>(0);
  constexpr char b = coef<p0, p1, p2>(1);
  constexpr char c = coef<p0, p1, p2>(2);
  constexpr char denom = b*b - a*c;
  if (i == 0) {
    return c/static_cast<double>(denom);
  } else if (i == 1) {
    return a/static_cast<double>(denom);
  } else {
    return -b/static_cast<double>(denom);
  }
}

template <char p0, char p1, char p2>
double eval_Q(double lam1, double lam2) {
  return coef<p0, p1, p2>(0)*lam1*lam1 + coef<p0, p1, p2>(1)*lam1*lam2 +
    coef<p0, p1, p2>(2)*lam2*lam2 + coef<p0, p1, p2>(3)*lam1 +
    coef<p0, p1, p2>(4)*lam2 + coef<p0, p1, p2>(5);
}

/**
 * Computes the descent direction for the Newton's method used for the
 * tetrahedral updates (i.e., the inverse of the Hessian times the
 * gradient of F0). This uses some simplifications to speed up
 * computations (see the paper).
 */
template <char p0, char p1, char p2>
std::pair<double, double> compute_g(double du1, double du2, double lam1,
                                    double lam2, double sh) {
  static constexpr char c0 = coef<p0, p1, p2>(0);
  static constexpr char c1 = coef<p0, p1, p2>(1);
  static constexpr char c2 = coef<p0, p1, p2>(2);
  static constexpr char c3 = coef<p0, p1, p2>(3);
  static constexpr char c4 = coef<p0, p1, p2>(4);
  double const alpha = c0*lam1 + c1*lam2 + c3;
  double const beta = c1*lam1 + c2*lam2 + c4;
  double const alpha_times_beta = alpha*beta;
  double const Q = eval_Q<p0, p1, p2>(lam1, lam2);
  double const A = c0 - alpha*alpha/Q;
  double const B = c1 - alpha_times_beta/Q;
  double const C = c2 - beta*beta/Q;
  double const scale = sqrt(Q)/(sh*(B*B - A*C));
  return {scale*(-C*du1 + B*du2), scale*(-B*du1 - A*du2)};
}

template <char p0, char p1, char p2>
static double compute_step_size(double du1, double du2, double lam1,
                                double lam2, double sh) {
  static constexpr double D[3] = {
    dPt_times_dP_coef<p0, p1, p2>(0),
    dPt_times_dP_coef<p0, p1, p2>(1),
    dPt_times_dP_coef<p0, p1, p2>(2)
  };
  static constexpr double q[2] = {
    D[0]*dpi_dot_p0<p0, p1, p2>(1) + D[1]*dpi_dot_p0<p0, p1, p2>(2),
    D[1]*dpi_dot_p0<p0, p1, p2>(1) + D[2]*dpi_dot_p0<p0, p1, p2>(2)
  };
  double const tmp1 = D[0]*du1*du1 + 2*D[1]*du1*du2 + D[2]*du2*du2;
  double const tmp2 = eval_Q<p0, p1, p2>(lam1, lam2);
  double const tmp3 = (lam1 + q[0])*du1 + (lam2 + q[1])*du2;
  double const beta = (tmp1 - tmp3*tmp3/tmp2)/(sh*sh);
  return -1/std::sqrt(1 - beta);
}

#endif // __UPDATE_RULES_TETRA_UPDATES_UTIL_HPP__
