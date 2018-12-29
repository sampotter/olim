#ifndef __UPDATES_TRI_IMPL_HPP__
#define __UPDATES_TRI_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#include <cmath>

#include "common.hpp"
#include "hybrid.hpp"
#include "updates.utils.hpp"

#define l__(x) std::sqrt((dp_dot_dp*(x) + 2*dp_dot_p0)*(x) + p0_dot_p0)
#define check__(x) std::fabs(alpha*l__(x) - dp_dot_p0 - dp_dot_dp*(x))

#define int_sqrt__(c) _sqrt_table[static_cast<int>(c)]
#define s__(x) (s + (1 - x)*s0 + x*s1)/2

#define F0_line__(i) (u##i + h*(s + s##i)*int_sqrt__(p##i##_dot_p##i)/2)

template <int n, int p0, int p1>
updates::info<1>
updates::tri_bv<MP0, n, p0, p1>::operator()(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  using std::min;

  // TODO: replace with constexpr function
  static constexpr double _sqrt_table[4] = {
    0.0,
    1.0,
    1.4142135623730951,
    1.7320508075688772
  };

  constexpr int p0_dot_p0 = dot(p0, p0);
  constexpr int p0_dot_p1 = dot(p0, p1);
  constexpr int p1_dot_p1 = dot(p1, p1);
  constexpr int dp_dot_p0 = p0_dot_p1 - p0_dot_p0;
  constexpr int dp_dot_p0_sq = dp_dot_p0*dp_dot_p0;
  constexpr int dp_dot_dp = p1_dot_p1 - 2*p0_dot_p1 + p0_dot_p0;

  double const du = u1 - u0;
  double const alpha = -du/(s__(0.5)*h), alpha_sq = alpha*alpha;
  double const tmp = alpha_sq - dp_dot_dp;
  double const a = dp_dot_dp*tmp;
  double const b = dp_dot_p0*tmp;
  double const c = alpha_sq*p0_dot_p0 - dp_dot_p0_sq;
  double const disc = b*b - a*c;

  info<1> info;
  if (disc < 0 || a == 0) {
    double const tmp0 = F0_line__(0), tmp1 = F0_line__(1);
    if (tmp0 < tmp1) {
      info.value = tmp0;
      info.lambda[0] = 0;
    } else {
      info.value = tmp1;
      info.lambda[0] = 1;
    }
  } else {
    double const lhs = -b/a, rhs = sqrt(disc)/a;
    double const lam1 = lhs - rhs, lam2 = lhs + rhs;
    info.lambda[0] = check__(lam1) < check__(lam2) ? lam1 : lam2;
    info.value = info.lambda[0] < 0 || 1 < info.lambda[0] ?
      min(F0_line__(0), F0_line__(1)) :
      u0 + info.lambda[0]*du + h*s__(info.lambda[0])*l__(info.lambda[0]);
  }

  return info;
}

#undef F0_line__

template <int n>
updates::info<1>
updates::tri<MP0, n>::operator()(
  double const * p0, double const * p1,
  double u0, double u1, double s, double s0, double s1, double h,
  double const * p_fac, double s_fac) const
{
  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  double sh = ((s + (s0 + s1)/2)/2)*h, shfac = s_fac*h;
  double lfac0 = dist2<n>(p0, p_fac);
  double lfac1 = dist2<n>(p1, p_fac);
  double T0 = shfac*lfac0, T1 = shfac*lfac1;
  double tau0 = u0 - T0, tau1 = u1 - T1, dtau = tau1 - tau0;

  double dp[n];
  sub<n>(p1, p0, dp);

  double nu[n], nufac[n];

  auto const grad = [&] (double lam) {
    axpy<n>(lam, dp, p0, nu);
    sub<n>(nu, p_fac, nufac);
    scal_inplace<n>(1./norm2<n>(nu), nu);
    scal_inplace<n>(1./norm2<n>(nufac), nufac);
    return dtau + shfac*dot<n>(dp, nufac) + sh*dot<n>(dp, nu);
  };

  double arglam;
  hybrid_status status;
  std::tie(arglam, status) = hybrid(grad, 0., 1.);

  info<1> info;
  if (status == hybrid_status::DEGENERATE) {
    double F0 = tau0 + T0 + sh*norm2<n>(p0);
    double F1 = tau1 + T1 + sh*norm2<n>(p1);
    info.lambda[0] = F0 < F1 ? 0 : 1;
    info.value = std::min(F0, F1);
  }
  else {
    info.lambda[0] = arglam;
    axpy<n>(arglam, dp, p0, nu);
    sub<n>(nu, p_fac, nufac);
    info.value = tau0 + dtau*arglam + shfac*norm2<n>(nufac) +
      sh*norm2<n>(nu);
  }
  return info;
}

template <int n>
updates::info<1>
updates::tri<MP0, n>::operator()(
  double const * p0, double const * p1, double u0, double u1,
  double s, double s0, double s1, double h) const
{
  using std::min;

  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  double dp[n], dp_dot_dp, dp_dot_p0, p0_dot_p0;
  sub<n>(p1, p0, dp);
  dp_dot_dp = dot<n>(dp, dp);
  dp_dot_p0 = dot<n>(dp, p0);
  p0_dot_p0 = dot<n>(p0, p0);

  double const du = u1 - u0;
  double const stheta = (s + (s0 + s1)/2)/2;
  double const alpha = -du/(stheta*h), alpha_sq = alpha*alpha;
  double const tmp = alpha_sq - dp_dot_dp;
  double const a = dp_dot_dp*tmp;
  double const b = dp_dot_p0*tmp;
  double const c = alpha_sq*p0_dot_p0 - dp_dot_p0*dp_dot_p0;
  double const disc = b*b - a*c;

  double const F0 = u0 + h*(s + s0)*norm2<n>(p0)/2;
  double const F1 = u1 + h*(s + s1)*norm2<n>(p1)/2;

  info<1> info;
  if (disc < 0 || a == 0) {
    if (F0 < F1) {
      info.value = F0;
      info.lambda[0] = 0;
    } else {
      info.value = F1;
      info.lambda[0] = 1;
    }
  } else {
    double const lhs = -b/a, rhs = sqrt(disc)/a;
    double const lam1 = lhs - rhs, lam2 = lhs + rhs;
    info.lambda[0] = check__(lam1) < check__(lam2) ? lam1 : lam2;
    info.value = info.lambda[0] < 0 || 1 < info.lambda[0] ?
      min(F0, F1) :
      u0 + info.lambda[0]*du + h*s__(info.lambda[0])*l__(info.lambda[0]);
  }
  return info;
}

#define F0_line__(i) (u##i + sh*int_sqrt__(p##i##_dot_p##i))

template <int n, int p0, int p1>
updates::info<1>
updates::tri_bv<RHR, n, p0, p1>::operator()(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  using std::min;

  (void) s0;
  (void) s1;

  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  static constexpr double _sqrt_table[4] = {
    0.0,
    1.0,
    1.4142135623730951,
    1.7320508075688772
  };

  constexpr int p0_dot_p0 = dot(p0, p0);
  constexpr int p0_dot_p1 = dot(p0, p1);
  constexpr int p1_dot_p1 = dot(p1, p1);
  constexpr int dp_dot_p0 = p0_dot_p1 - p0_dot_p0;
  constexpr int dp_dot_p0_sq = dp_dot_p0*dp_dot_p0;
  constexpr int dp_dot_dp = p1_dot_p1 - 2*p0_dot_p1 + p0_dot_p0;

  double const du = u1 - u0;
  double const sh = s*h;
  double const alpha = -du/sh, alpha_sq = alpha*alpha;
  double const tmp = alpha_sq - dp_dot_dp;
  double const a = dp_dot_dp*tmp;
  double const b = dp_dot_p0*tmp;
  double const c = alpha_sq*p0_dot_p0 - dp_dot_p0_sq;
  double const disc = b*b - a*c;

  info<1> info;
  if (disc < 0 || a == 0) {
    double const tmp0 = F0_line__(0), tmp1 = F0_line__(1);
    if (tmp0 < tmp1) {
      info.value = tmp0;
      info.lambda[0] = 0;
    } else {
      info.value = tmp1;
      info.lambda[0] = 1;
    }
  } else {
    double const lhs = -b/a, rhs = sqrt(disc)/a;
    double const lam1 = lhs - rhs, lam2 = lhs + rhs;
    info.lambda[0] = check__(lam1) < check__(lam2) ? lam1 : lam2;
    info.value = info.lambda[0] < 0 || 1 < info.lambda[0] ?
      min(F0_line__(0), F0_line__(1)) :
      u0 + info.lambda[0]*du + sh*l__(info.lambda[0]);
  }
  return info;
}

#undef int_sqrt__
#undef s__

template <int n>
updates::info<1>
updates::tri<RHR, n>::operator()(
  double const * p0, double const * p1,
  double u0, double u1, double s, double s0, double s1, double h,
  double const * p_fac, double s_fac) const
{
  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  (void) s0;
  (void) s1;

  double sh = s*h, shfac = s_fac*h;
  double lfac0 = dist2<n>(p0, p_fac);
  double lfac1 = dist2<n>(p1, p_fac);
  double T0 = shfac*lfac0, T1 = shfac*lfac1;
  double tau0 = u0 - T0, tau1 = u1 - T1, dtau = tau1 - tau0;

  double dp[n];
  sub<n>(p1, p0, dp);

  double nu[n], nufac[n];

  auto const grad = [&] (double lam) {
    axpy<n>(lam, dp, p0, nu);
    sub<n>(nu, p_fac, nufac);
    scal_inplace<n>(1./norm2<n>(nu), nu);
    scal_inplace<n>(1./norm2<n>(nufac), nufac);
    return dtau + shfac*dot<n>(dp, nufac) + sh*dot<n>(dp, nu);
  };

  double arglam;
  hybrid_status status;
  std::tie(arglam, status) = hybrid(grad, 0., 1.);

  info<1> info;
  if (status == hybrid_status::DEGENERATE) {
    double F0 = tau0 + T0 + sh*norm2<n>(p0);
    double F1 = tau1 + T1 + sh*norm2<n>(p1);
    info.lambda[0] = F0 < F1 ? 0 : 1;
    info.value = std::min(F0, F1);
  }
  else {
    info.lambda[0] = arglam;
    axpy<n>(arglam, dp, p0, nu);
    sub<n>(nu, p_fac, nufac);
    info.value = tau0 + dtau*arglam + shfac*norm2<n>(nufac) +
      sh*norm2<n>(nu);
  }
  return info;
}

template <int n>
updates::info<1>
updates::tri<RHR, n>::operator()(
  double const * p0, double const * p1,
  double u0, double u1, double s, double s0, double s1, double h) const
{
  (void) s0;
  (void) s1;

  using std::min;

  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  double dp[n], dp_dot_dp, dp_dot_p0, p0_dot_p0;
  sub<n>(p1, p0, dp);
  dp_dot_dp = dot<n>(dp, dp);
  dp_dot_p0 = dot<n>(dp, p0);
  p0_dot_p0 = dot<n>(p0, p0);

  double const du = u1 - u0;
  double const sh = s*h;
  double const alpha = -du/sh, alpha_sq = alpha*alpha;
  double const tmp = alpha_sq - dp_dot_dp;
  double const a = dp_dot_dp*tmp;
  double const b = dp_dot_p0*tmp;
  double const c = alpha_sq*p0_dot_p0 - dp_dot_p0*dp_dot_p0;
  double const disc = b*b - a*c;

  double const F0 = u0 + sh*norm2<n>(p0);
  double const F1 = u1 + sh*norm2<n>(p1);

  info<1> info;
  if (disc < 0 || a == 0) {
    if (F0 < F1) {
      info.value = F0;
      info.lambda[0] = 0;
    } else {
      info.value = F1;
      info.lambda[0] = 1;
    }
  } else {
    double const lhs = -b/a, rhs = sqrt(disc)/a;
    double const lam1 = lhs - rhs, lam2 = lhs + rhs;
    info.lambda[0] = check__(lam1) < check__(lam2) ? lam1 : lam2;
    info.value = info.lambda[0] < 0 || 1 < info.lambda[0] ?
      min(F0, F1) :
      u0 + info.lambda[0]*du + sh*l__(info.lambda[0]);
  }
  return info;
}

#undef check__
#undef l__

#define u__(x) ((1 - (x))*u0 + (x)*u1)
#define q__(x) ((dp_dot_dp*x + 2*dp_dot_p0)*(x) + p0_dot_p0)
#define l__(x) std::sqrt(q__(x))
#define s__(x) ((s + (1 - (x))*s0 + (x)*s1)/2)
#define F1__(x) (u__(x) + h*s__(x)*l__(x))
#define dF1__(x) (du + h*(ds*q__(x)/2 + s__(x)*dp_dot_plam)/l__(x))
#define d2F1__(x) h*(s__(x)*dp_dot_plam*dp_dot_plam/q__(x) + \
                     ds*dp_dot_plam/2 + 2*s__(x)*dp_dot_dp)/l__(x)

/**
 * F1 specialization
 */
template <int n, int p0, int p1>
updates::info<1>
updates::tri_bv<MP1, n, p0, p1>::operator()(
  double u0, double u1, double s, double s0, double s1, double h) const
{
  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  constexpr int p0_dot_p0 = dot(p0, p0);
  constexpr int p0_dot_p1 = dot(p0, p1);
  constexpr int p1_dot_p1 = dot(p1, p1);
  constexpr int dp_dot_p0 = p0_dot_p1 - p0_dot_p0;
  constexpr int dp_dot_dp = p1_dot_p1 - 2*p0_dot_p1 + p0_dot_p0;

  double const ds = s1 - s0, du = u1 - u0;
  double dp_dot_p_lam, l_lam, s_lam;

  auto const grad = [&] (double lam) {
    dp_dot_p_lam = dp_dot_p0 + lam*dp_dot_dp;
    l_lam = std::sqrt(p0_dot_p0 + lam*(dp_dot_p0 + dp_dot_p_lam));
    s_lam = (s + s0 + ds*lam)/2;
    return du + h*(l_lam*ds/2 + s_lam*dp_dot_p_lam/l_lam);
  };

  double arglam;
  hybrid_status status;
  std::tie(arglam, status) = hybrid(grad, 0., 1.);
  
  info<1> info;
  if (status == hybrid_status::DEGENERATE) {
    double F0 = u0 + (s + s0)*h*std::sqrt(p0_dot_p0)/2;
    double F1 = u1 + (s + s1)*h*std::sqrt(p1_dot_p1)/2;
    info.lambda[0] = F0 < F1 ? 0 : 1;
    info.value = std::min(F0, F1);
  }
  else {
    info.lambda[0] = arglam;
    dp_dot_p_lam = dp_dot_p0 + arglam*dp_dot_dp;
    l_lam = std::sqrt(p0_dot_p0 + arglam*(dp_dot_p0 + dp_dot_p_lam));
    s_lam = (s + s0 + ds*arglam)/2;
    info.value = u0 + du*arglam + s_lam*h*l_lam;
  }
  return info;
}

template <int n>
updates::info<1>
updates::tri<MP1, n>::operator()(
  double const * p0, double const * p1,
  double u0, double u1, double s, double s0, double s1, double h,
  double const * p_fac, double s_fac) const
{
  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  double shfac = s_fac*h;
  double ds = s1 - s0;
  double T0 = shfac*dist2<n>(p0, p_fac);
  double T1 = shfac*dist2<n>(p1, p_fac);
  double tau0 = u0 - T0, tau1 = u1 - T1, dtau = tau1 - tau0;

  double dp[n], nu[n], nufac[n];
  sub<n>(p1, p0, dp);

  double s_lam, l_lam;

  auto const grad = [&] (double lam) {
    axpy<n>(lam, dp, p0, nu);
    sub<n>(nu, p_fac, nufac);
    l_lam = norm2<n>(nu);
    scal_inplace<n>(1./l_lam, nu);
    scal_inplace<n>(1./norm2<n>(nufac), nufac);
    s_lam = (s + s0 + ds*lam)/2;
    return dtau + shfac*dot<n>(dp, nufac) +
      h*(s_lam*dot<n>(dp, nu) + l_lam*ds/2);
  };

  double arglam;
  hybrid_status status;
  std::tie(arglam, status) = hybrid(grad, 0., 1.);

  info<1> info;
  if (status == hybrid_status::DEGENERATE) {
    double F0 = tau0 + T0 + (s + s0)*h*norm2<n>(p0)/2;
    double F1 = tau1 + T1 + (s + s1)*h*norm2<n>(p1)/2;
    info.lambda[0] = F0 < F1 ? 0 : 1;
    info.value = std::min(F0, F1);
  }
  else {
    info.lambda[0] = arglam;
    axpy<n>(arglam, dp, p0, nu);
    sub<n>(nu, p_fac, nufac);
    s_lam = (s + s0 + ds*arglam)/2;
    info.value = tau0 + dtau*arglam + shfac*norm2<n>(nufac) +
      s_lam*h*norm2<n>(nu);
  }
  return info;
}

template <int n>
updates::info<1>
updates::tri<MP1, n>::operator()(
  double const * p0, double const * p1, double u0, double u1,
  double s, double s0, double s1, double h) const
{
  assert(s > 0);
  assert(s0 > 0);
  assert(s1 > 0);

  double dp[n], dp_dot_dp, dp_dot_p0, p0_dot_p0;
  sub<n>(p1, p0, dp);
  dp_dot_dp = dot<n>(dp, dp);
  dp_dot_p0 = dot<n>(dp, p0);
  p0_dot_p0 = dot<n>(p0, p0);

  constexpr double c1 = 1e-4;

  double const ds = s1 - s0, du = u1 - u0;

  info<1> info;

  // Check gradients at endpoints to see if we can skip this update
  double dp_dot_plam = dp_dot_p0;
  if (dF1__(0) > 0) {
    info.value = F1__(0);
    info.lambda[0] = 0;
    return info;
  }
  dp_dot_plam += dp_dot_dp;
  if (dF1__(1) < 0) {
    info.value = F1__(1);
    info.lambda[0] = 1;
    return info;
  }

  // Try to minimize F1 by starting from the mp0 minimizer using
  // Newton's method. This may fail if the function isn't
  // well-behaved. To try to determine when this happens, keep track
  // of the iteration count, and use a modified Newton's method if we
  // exceed a limit of 10 iterations.
  // 
  // TODO: this is messy---should replace this with a Newton iteration
  // in a separate file
  info = tri<MP0, n>()(p0, p1, u0, u1, s, s0, s1, h);
  double lam = info.lambda[0], g, iter = 0;
  do {
    dp_dot_plam = dp_dot_p0 + lam*dp_dot_dp;
    g = -dF1__(lam)/d2F1__(lam);
    lam = std::max(0., std::min(1., lam + g));
  } while (iter++ < 10 && fabs(g) > eps<double>);

  if (iter == 10) {
    bool conv;
    double lam[2], F1[2], dF1, d2F1, g, alpha;
    lam[0] = info.lambda[0];
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
      conv = fabs(lam[1] - lam[0]) <= eps<double> ||
        fabs(F1[1] - F1[0]) <= eps<double>;
      lam[0] = lam[1];
      F1[0] = F1[1];
    } while (!conv);
    info.value = F1[1];
    info.lambda[0] = lam[1];
  } else {
    info.value = F1__(lam);
    info.lambda[0] = lam;
  }
  return info;
}

#undef u__
#undef q__
#undef l__
#undef theta__
#undef s__
#undef F1__
#undef dF1__
#undef d2F1__

#endif // __UPDATES_TRI_IMPL_HPP__
