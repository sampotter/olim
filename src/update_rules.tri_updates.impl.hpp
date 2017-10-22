#ifndef __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__

#include <src/config.hpp>

#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.defs.hpp"
#include "common.macros.hpp"
#ifdef EIKONAL_DEBUG
#    include "olim_util.hpp"
#endif

namespace update_rules {
  template <class speed_estimator, bool is_constrained>
  std::enable_if_t<!is_constrained, double>
  rect_tri_updates<speed_estimator, is_constrained>::tri11(
    double u0, double u1, double s, double s0, double s1, double h) const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, h, s, s0, s1);
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
    printf("tri11(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
    return T;
  }

  template <class speed_estimator, bool is_constrained>
  std::enable_if_t<!is_constrained, double>
  rect_tri_updates<speed_estimator, is_constrained>::tri12(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, h, s, s0, s1);
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
    printf("tri12(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
    return T;
  }

  template <class speed_estimator, bool is_constrained>
  std::enable_if_t<!is_constrained, double>
  rect_tri_updates<speed_estimator, is_constrained>::tri13(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, h, s, s0, s1);
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
    printf("tri13(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
    return T;
  }

  template <class speed_estimator, bool is_constrained>
  std::enable_if_t<!is_constrained, double>
  rect_tri_updates<speed_estimator, is_constrained>::tri22(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, h, s, s0, s1);
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
    printf("tri22(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
    return T;
  }

  template <class speed_estimator, bool is_constrained>
  std::enable_if_t<!is_constrained, double>
  rect_tri_updates<speed_estimator, is_constrained>::tri23(
    double u0, double u1, double s, double s0, double s1, double h)
    const
  {
#ifdef EIKONAL_DEBUG
    check_params(u0, u1, h, s, s0, s1);
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
    printf("tri23(u0 = %g, u1 = %g, s = %g, h = %g) -> %g\n", u0, u1, s, h, T);
#endif
    return T;
  }
}

#endif // __UPDATE_RULES_TRI_UPDATES_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
