#ifndef __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__

#include <src/config.hpp>

#include <random>

#if EIKONAL_DEBUG
#    include <cassert>
#endif

#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.hpp"
#include "common.defs.hpp"
#include "common.macros.hpp"
#include "olim_util.hpp"
#include "qroots.hpp"

namespace update_rules {
  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri11(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
    return tri11_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri12(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
    return tri12_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri13(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
    return tri13_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri22(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
    return tri22_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri23(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
    return tri23_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri11_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::true_type &&) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double sh = this->s_hat(s, s0, s1)*h;
    double c = (u0 - u1)/sh;
    double lam = 0.5 + (c > 0 ? 1 : -1)*std::fabs(c)/(2*sqrt(2 - c*c));
#if PRINT_UPDATES
    double tmp = (1 - lam)*u0 + lam*u1 + sh*sqrt(2*lam*(lam - 1) + 1);
    printf("tri11(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g, "
           "is_constrained = true) -> %g\n",
           u0, u1, s, s0, s1, h, tmp);
    return tmp;
#endif
    return (1 - lam)*u0 + lam*u1 + s*h*sqrt(2*lam*(lam - 1) + 1);
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri11_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double sh = this->s_hat(s, s0, s1)*h, du = u1 - u0, alpha = fabs(du)/sh,
      alpha_sq = alpha*alpha;
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
    printf("tri11(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g, "
           "is_constrained = false) -> %g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri12_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::true_type &&) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double sh = this->s_hat(s, s0, s1)*h, c = std::fabs(u0 - u1)/sh;
    double sgn = u0 > u1 ? 1 : -1;
    double lam = std::max(0.0, std::min(1.0, sgn*c/sqrt(1 - c*c)));
#if PRINT_UPDATES
    double tmp = (1 - lam)*u0 + lam*u1 + sh*std::sqrt(lam*lam + 1);
    printf("tri12(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g, "
           "is_constrained = true) -> %g\n", u0, u1, s, s0, s1, h, tmp);
    return tmp;
#endif
    return (1 - lam)*u0 + lam*u1 + sh*std::sqrt(lam*lam + 1);
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri12_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&)
    const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double sh = this->s_hat(s, s0, s1)*h, du = u1 - u0, alpha = fabs(du)/sh,
      alpha_sq = alpha*alpha;
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
    printf("tri12(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g, "
           "is_constrained = false) -> %g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri13_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&)
    const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double sh = this->s_hat(s, s0, s1)*h, du = u1 - u0, alpha = fabs(du)/sh,
      alpha_sq = alpha*alpha;
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
    printf("tri13(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g, "
           "is_constrained = false) -> %g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri22_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&)
    const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double sh = this->s_hat(s, s0, s1)*h, du = u1 - u0, alpha = fabs(du)/sh,
      alpha_sq = alpha*alpha;
    double T = INF(double);
    if (alpha_sq > 2) {
      return T;
    }
    double rhs = sqrt3*alpha/sqrt(2 - alpha_sq);
    double lam = (1 + rhs)/2, l = sqrt(2*(1 - lam + lam*lam));
    if (0 <= lam && lam <= 1 && fabs(du/sh + (2*lam - 1)/l) < 1e-13) {
      T = (1 - lam)*u0 + lam*u1 + sh*l;
    } else {
      lam = (1 - rhs)/2, l = sqrt(2*(1 - lam + lam*lam));
      if (0 <= lam && lam <= 1 && fabs(du/sh + (2*lam - 1)/l) < 1e-13) {
        T = (1 - lam)*u0 + lam*u1 + sh*l;
      }
    }
#if PRINT_UPDATES
    printf("tri22(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g, "
           "is_constrained = false) -> %g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

  template <class speed_estimator, bool is_constrained>
  double
  rect_tri_updates<speed_estimator, is_constrained>::tri23_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&)
    const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double sh = this->s_hat(s, s0, s1)*h, du = u1 - u0, alpha = fabs(du)/sh,
      alpha_sq = alpha*alpha;
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
    printf("tri23(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g, "
           "is_constrained = false) -> %g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

#define EVAL_q(lam) ((A*(lam) + 2*B)*(lam) + C)
#define EVAL_l(lam) std::sqrt(EVAL_q(lam))
#define EVAL_sbar(lam) (sbar0 + (lam)*dsbar)
#define EVAL_u(lam) (u0 + (lam)*du)

  template <int A, int B, int C>
  inline double mp1_tri_newton(double u0, double u1, double s, double s0,
                               double s1, double h) {
    using std::max;
    using std::min;

    double const sbar0 = (s + s0)/2;
    double const sbar1 = (s + s1)/2;
    double const dsbar = sbar1 - sbar0;
    double const du = u1 - u0;

    auto const d2F = [&] (double lam) -> double {
      double q = EVAL_q(lam);
      double l = std::sqrt(q);
      double dl = (A*lam + B)/l;
      double d2l = (A*C - B*B)/(q*l);
      return h*(2*dsbar*dl + EVAL_sbar(lam)*d2l);
    };

    static std::random_device dev;
    static std::minstd_rand0 gen(dev());
    static std::uniform_real_distribution<double> dist(0.4, 0.6);

    bool found_minima = false;
    double lam, F0, F1, dF, dlam, alpha, g, l, u, sbar, c1 = 1e-4, eps = 10*EPS(double);
    bool out_of_bounds = false;

    while (!found_minima) {
      lam = dist(gen);

      // Initialize values for iteration
      l = EVAL_l(lam);
      u = EVAL_u(lam);
      sbar = EVAL_sbar(lam);
      F1 = u + h*sbar*l;
      do {
        F0 = F1; // Save old F value

        // Compute the first derivative of F
        dF = du + h*(dsbar*l + sbar*(A*lam + B)/l);

        // Compute descent direction
        dlam = -dF/max(0.1, d2F(lam));

        // Backtracking line search to satisfy Wolfe conditions
        g = c1*dF*dlam;
        alpha = 1;
        double tmp = lam + alpha*dlam;
        double lhs = EVAL_u(tmp) + h*EVAL_sbar(tmp)*EVAL_l(tmp);
        while (lhs - eps > F1 + alpha*g) {
          alpha *= 0.99;
          tmp = lam + alpha*dlam;
          lhs = EVAL_u(tmp) + h*EVAL_sbar(tmp)*EVAL_l(tmp);
        }

        // Update lambda iterate
        lam += alpha*dlam;

        // Keep track of whether the iterate has exited the feasible
        // set (gone `out of bounds'). If it's been `out of bounds'
        // for for two iterations in a row, return infinity (we're
        // solving the constrained problem)
        if (lam < 0 || lam >= 1) {
          if (out_of_bounds) {
            return INF(double);
          }
          out_of_bounds = true;
        }

        // Update iteration's values
        l = EVAL_l(lam);
        u = EVAL_u(lam);
        sbar = EVAL_sbar(lam);
        F1 = u + h*sbar*l;
      } while (fabs(dlam) > eps*fabs(dlam) && fabs(F1 - F0) > eps);
      if (lam < 0 || lam >= 1) {
        return INF(double);
      }
      if (d2F(lam) > 0) {
        found_minima = true;
      }
    }
    return F1;
  }

  template <bool is_constrained>
  double mp1_tri_updates<is_constrained>::tri11(
    double u0, double u1, double s, double s0, double s1, double h) const
  {
    return tri11_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <bool is_constrained>
  double mp1_tri_updates<is_constrained>::tri12(
    double u0, double u1, double s, double s0, double s1, double h) const
  {
    return tri12_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <bool is_constrained>
  double mp1_tri_updates<is_constrained>::tri13(
    double u0, double u1, double s, double s0, double s1, double h) const
  {
    return tri13_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <bool is_constrained>
  double mp1_tri_updates<is_constrained>::tri22(
    double u0, double u1, double s, double s0, double s1, double h) const
  {
    return tri22_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

  template <bool is_constrained>
  double mp1_tri_updates<is_constrained>::tri23(
    double u0, double u1, double s, double s0, double s1, double h) const
  {
    return tri23_impl(
      u0, u1, s, s0, s1, h, eikonal::bool_t<is_constrained> {});
  }

#ifdef EIKONAL_DEBUG
#    if !RELWITHDEBINFO
  template <int A, int B, int C>
  static inline void
  verify_edge_solution(double u0, double u1, double s, double s0, double s1,
                       double h) {
    using std::sqrt;
    double const sbar0 = (s + s0)/2, sbar1 = (s + s1)/2;
    double const l0 = sqrt(C), l1 = sqrt(A + 2*B + C);
    double const T_left = u0 + h*sbar0*l0, T_right = u1 + h*sbar1*l1;
    double const dlam = sqrt(EPS(double));
    double const lam1 = T_left < T_right ? dlam : 1 - dlam;
    double const lam0 = 1 - lam1;
    double const l = sqrt(A*lam1*lam1 + 2*B*lam1 + C);
    double const T = lam0*u0 + lam1*u1 + h*(sbar0*lam0 + sbar1*lam1)*l;
    (void) T;
    if (T_left < T_right) {
      assert(T_left <= T);
    } else {
      assert(T_right <= T);
    }
  }
#    endif
#endif

  template <bool is_constrained>
  double
  mp1_tri_updates<is_constrained>::tri11_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::true_type &&) const
  {
    double T = tri11_impl(u0, u1, s, s0, s1, h, std::false_type {});
#ifdef EIKONAL_DEBUG
#    if !RELWITHDEBINFO
    if (std::isinf(T)) {
      verify_edge_solution<2, -1, 1>(u0, u1, s, s0, s1, h);
    }
#    endif
#endif
    return std::isinf(T) ?
      std::min(u0 + (s + s0)*h/2, u1 + (s + s1)*h/2) : T;
  }

  template <bool is_constrained>
  double
  mp1_tri_updates<is_constrained>::tri11_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double T = mp1_tri_newton<2, -1, 1>(u0, u1, s, s0, s1, h);
#ifdef PRINT_UPDATES
    printf("tri11(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g) -> "
           "%g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

  template <bool is_constrained>
  double
  mp1_tri_updates<is_constrained>::tri12_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::true_type &&) const
  {
    double T = tri12_impl(u0, u1, s, s0, s1, h, std::false_type {});
#ifdef EIKONAL_DEBUG
#if !RELWITHDEBINFO
    if (std::isinf(T)) {
      verify_edge_solution<1, 0, 1>(u0, u1, s, s0, s1, h);
    }
#    endif
#endif
    return std::isinf(T) ?
      std::min(u0 + (s + s0)*h/2, u1 + (s + s1)*h*sqrt2/2) : T;
  }

  template <bool is_constrained>
  double
  mp1_tri_updates<is_constrained>::tri12_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double T = mp1_tri_newton<1, 0, 1>(u0, u1, s, s0, s1, h);
#if PRINT_UPDATES
    printf("tri12(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g) -> "
           "%g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

  template <bool is_constrained>
  double
  mp1_tri_updates<is_constrained>::tri13_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double T = mp1_tri_newton<2, 0, 1>(u0, u1, s, s0, s1, h);
#ifdef PRINT_UPDATES
    printf("tri13(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g) -> "
           "%g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

  template <bool is_constrained>
  double
  mp1_tri_updates<is_constrained>::tri22_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double T = mp1_tri_newton<2, -1, 2>(u0, u1, s, s0, s1, h);
#ifdef PRINT_UPDATES
    printf("tri22(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g) -> "
           "%g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }

  template <bool is_constrained>
  double
  mp1_tri_updates<is_constrained>::tri23_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::false_type &&) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, s, s0, s1, h);
#endif
    double T = mp1_tri_newton<1, 0, 2>(u0, u1, s, s0, s1, h);
#ifdef PRINT_UPDATES
    printf("tri23(u0 = %g, u1 = %g, s = %g, s0 = %g, s1 = %g, h = %g) -> "
           "%g\n", u0, u1, s, s0, s1, h, T);
#endif
    return T;
  }
}

#endif // __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__
