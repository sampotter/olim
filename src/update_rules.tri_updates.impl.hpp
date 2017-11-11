#ifndef __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__

#include <src/config.hpp>

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

  template <int A, int B, int C>
  inline double mp1_tri_newton(double u0, double u1, double s, double s0,
                               double s1, double h) {
    double const sbar0 = (s + s0)/2;
    double const sbar1 = (s + s1)/2;
    double const dsbar = sbar1 - sbar0;
    double const du = u1 - u0;

    auto const u = [&] (double lam) -> double {
      return (1 - lam)*u0 + lam*u1;
    };

    auto const sbar = [&] (double lam) -> double {
      return (1 - lam)*sbar0 + lam*sbar1;
    };

    auto const q = [&] (double lam) -> double {
      return A*lam*lam + 2*B*lam + C;
    };

    auto const dq = [&] (double lam) -> double {
      return 2*(A*lam + B);
    };

    auto const l = [&] (double lam) -> double {
      return std::sqrt(q(lam));
    };

    auto const dl = [&] (double lam) -> double {
      return dq(lam)/(2*l(lam));
    };

    auto const d2l = [&] (double lam) -> double {
      return (A*C - B*B)/(q(lam)*l(lam));
    };

    auto const F = [&] (double lam) -> double {
      return u(lam) + h*sbar(lam)*l(lam);
    };

    auto const dF = [&] (double lam) -> double {
      return du + h*(dsbar*l(lam) + sbar(lam)*dl(lam));
    };

    auto const d2F = [&] (double lam) -> double {
      return h*(2*dsbar*dl(lam) + sbar(lam)*d2l(lam));
    };

    auto const p = [&] (double lam) -> double {
      return dF(lam)/d2F(lam);
    };

    double lam = 0.5, F0, F1 = F(lam), dlam;
    do {
      F0 = F1;
      dlam = p(lam);
      lam -= dlam;
      if (lam < 0 || 1 < lam) {
        return INF(double);
      }
      F1 = F(lam);
    } while (fabs(dlam)/fabs(lam) > EPS(double) && fabs(F1 - F0) > EPS(double));

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
    if (T_left < T_right) {
      assert(T_left <= T);
    } else {
      assert(T_right <= T);
    }
  }
#endif

  template <bool is_constrained>
  double
  mp1_tri_updates<is_constrained>::tri11_impl(
    double u0, double u1, double s, double s0, double s1, double h,
    std::true_type &&) const
  {
    double T = tri11_impl(u0, u1, s, s0, s1, h, std::false_type {});
#ifdef EIKONAL_DEBUG
    if (std::isinf(T)) {
      verify_edge_solution<2, -1, 1>(u0, u1, s, s0, s1, h);
    }
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
    if (std::isinf(T)) {
      verify_edge_solution<1, 0, 1>(u0, u1, s, s0, s1, h);
    }
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
