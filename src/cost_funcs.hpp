#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

#include <cassert>
#include <cmath>
#include <type_traits>

#include "bitops.hpp"
#include "common.macros.hpp"
#include "updates.utils.hpp"

enum class cost_func {mp0, mp1, rhr};

constexpr auto MP0 = cost_func::mp0;
constexpr auto MP1 = cost_func::mp1;
constexpr auto RHR = cost_func::rhr;

constexpr int sym_mat_size(int d) {
  return ((d + 1)*d)/2;
}

template <int d>
void lagmults(double const * lam, double const * df, double const * d2f,
              double * mu, int * k);

template <int d>
struct eval_wkspc {
  double u0, du[d], u_lam, sh_lam, l_lam;
};

template <int d>
struct fac_wkspc {
  double tau0, dtau[d], tau_lam, sh_lam, l_lam, sh_fac, l_fac_lam,
    dPt_nu_fac_lam[d];
};

template <int d, class base_wkspc = eval_wkspc<d>>
struct F0_wkspc: public base_wkspc {
  double dPt_nu_lam[d], dPt_dP[sym_mat_size(d)];
};

template <int d, class base_wkspc = eval_wkspc<d>>
struct F1_wkspc: public F0_wkspc<d, base_wkspc> {
  double sh_bar, theta_h_ds[d];
};

template <cost_func F, int d>
using F_wkspc = std::conditional_t<F == MP1, F1_wkspc<d>, F0_wkspc<d>>;

template <int d>
using F0_fac_wkspc = F0_wkspc<d, fac_wkspc<d>>;

template <int d>
using F1_fac_wkspc = F1_wkspc<d, fac_wkspc<d>>;

template <cost_func F, int d>
using F_fac_wkspc = std::conditional_t<
  F == MP1, F1_fac_wkspc<d>, F0_fac_wkspc<d>>;

template <cost_func F, class wkspc>
void set_sh_lam(wkspc & w, double s, double s0, double s1, double s2, double h)
{
  if (F == MP0) {
    w.sh_lam = h*(s + (s0 + s1 + s2)/3)/2;
  } else if (F == MP1 || F == RHR) {
    // NOTE: MP1 and RHR do different things here, and treat the
    // variable sh_lam differently. Later, MP1 also computes a
    // perturbation which is used when eval is called, and treats
    // sh_lam as the base term.
    w.sh_lam = h*s;
  }
}

template <cost_func F>
void set_args_common(F0_wkspc<2> & w, double u0, double u1, double u2,
                     double s, double s0, double s1, double s2, double h)
{
  w.u0 = u0;
  w.du[0] = u1 - u0;
  w.du[1] = u2 - u0;

  set_sh_lam<F, decltype(w)>(w, s, s0, s1, s2, h);
}

template <cost_func F>
void set_args_common(F1_wkspc<2> & w, double u0, double u1, double u2,
                     double s, double s0, double s1, double s2, double h)
{
  static_assert(F == MP1, "Cost function must be MP1");

  set_args_common<F>(
    static_cast<F0_wkspc<2> &>(w), u0, u1, u2, s, s0, s1, s2, h);

  w.sh_bar = (s + s0)*h/2;

  w.theta_h_ds[0] = h*(s1 - s0)/2;
  w.theta_h_ds[1] = h*(s2 - s0)/2;
}

template <cost_func F, int n, int p0, int p1, int p2>
void set_args(F_wkspc<F, 2> & w, double u0, double u1, double u2,
              double s, double s0, double s1, double s2, double h)
{
  using namespace bitops;
  using dim_t = dim<n>;

  set_args_common<F>(w, u0, u1, u2, s, s0, s1, s2, h);

  w.dPt_dP[0] = p_dot_q<p1, p1>(dim_t {}) - 2*p_dot_q<p1, p0>(dim_t {})
    + p_dot_q<p0, p0>(dim_t {});
  w.dPt_dP[1] = p_dot_q<p1, p2>(dim_t {}) - p_dot_q<p1, p0>(dim_t {})
    - p_dot_q<p2, p0>(dim_t {}) + p_dot_q<p0, p0>(dim_t {});
  w.dPt_dP[2] = p_dot_q<p2, p2>(dim_t {}) - 2*p_dot_q<p2, p0>(dim_t {})
    + p_dot_q<p0, p0>(dim_t {});
}

template <int n, class wkspc>
void set_dPt_dP(wkspc & w, double const * p0, double const * p1,
                double const * p2)
{
  double dp1[n], dp2[n];
  sub<n>(p1, p0, dp1);
  sub<n>(p2, p0, dp2);

  w.dPt_dP[0] = dot<n>(dp1, dp1);
  w.dPt_dP[1] = dot<n>(dp1, dp2);
  w.dPt_dP[2] = dot<n>(dp2, dp2);
}

template <cost_func F, int n>
void set_args(F_wkspc<F, 2> & w,
              double const * p0, double const * p1, double const * p2,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h)
{
  set_args_common<F>(w, u0, u1, u2, s, s0, s1, s2, h);
  set_dPt_dP<n>(w, p0, p1, p2);
}

template <cost_func F, int n>
void set_args(F0_wkspc<2, fac_wkspc<2>> & w,
              double const * p0, double const * p1, double const * p2,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h,
              double const * p_fac, double s_fac)
{
  w.sh_fac = s_fac*h;

  double const T0 = w.sh_fac*dist2<n>(p0, p_fac);

  w.tau0 = u0 - T0;
  w.dtau[0] = u1 - u0 - (w.sh_fac*dist2<n>(p1, p_fac) - T0);
  w.dtau[1] = u2 - u0 - (w.sh_fac*dist2<n>(p2, p_fac) - T0);

  set_sh_lam<F>(w, s, s0, s1, s2, h);

  set_dPt_dP<n>(w, p0, p1, p2);
}

template <cost_func F, int n>
void set_args(F1_wkspc<2, fac_wkspc<2>> & w,
              double const * p0, double const * p1, double const * p2,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h,
              double const * p_fac, double s_fac)
{
  static_assert(F == MP1, "Cost function must be MP1");

  set_args<F, n>(
    static_cast<F0_wkspc<2, fac_wkspc<2>> &>(w),
    p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, p_fac, s_fac);

  w.sh_bar = (s + s0)*h/2;

  w.theta_h_ds[0] = h*(s1 - s0)/2;
  w.theta_h_ds[1] = h*(s2 - s0)/2;
}

inline void set_lambda_common(F0_wkspc<2> & w, double const * lam)
{
  w.u_lam = w.u0 + w.du[0]*lam[0] + w.du[1]*lam[1];
}

inline void set_lambda_common(F1_wkspc<2> & w, double const * lam)
{
  set_lambda_common(static_cast<F0_wkspc<2> &>(w), lam);

  w.sh_lam = w.sh_bar + w.theta_h_ds[0]*lam[0] + w.theta_h_ds[1]*lam[1];
}

template <cost_func F, int n, int p0, int p1, int p2>
void set_lambda(F_wkspc<F, 2> & w, double const * lam)
{
  using namespace bitops;
  using dim_t = dim<n>;

  set_lambda_common(w, lam);

  w.l_lam = std::sqrt(
    p_dot_q<p0, p0>(dim_t {}) +
    2*(p_dot_q<p1, p0>(dim_t {}) - p_dot_q<p0, p0>(dim_t {}))*lam[0] +
    2*(p_dot_q<p2, p0>(dim_t {}) - p_dot_q<p0, p0>(dim_t {}))*lam[1] +
    dPt_dP<p0, p1, p2, 0>(dim_t {})*lam[0]*lam[0] +
    2*dPt_dP<p0, p1, p2, 1>(dim_t {})*lam[0]*lam[1] +
    dPt_dP<p0, p1, p2, 2>(dim_t {})*lam[1]*lam[1]
    // w.dPt_dP[0]*lam[0]*lam[0] +
    // 2*w.dPt_dP[1]*lam[0]*lam[1] +
    // w.dPt_dP[2]*lam[1]*lam[1]
    );
  CHECK(w.l_lam);

  w.dPt_nu_lam[0] = (
    p_dot_q<p1, p0>(dim_t {}) - p_dot_q<p0, p0>(dim_t {}) +
    dPt_dP<p0, p1, p2, 0>(dim_t {})*lam[0] +
    dPt_dP<p0, p1, p2, 1>(dim_t {})*lam[1]
    // w.dPt_dP[0]*lam[0] +
    // w.dPt_dP[1]*lam[1]
  )/w.l_lam;
  CHECK(w.dPt_nu_lam[0]);

  w.dPt_nu_lam[1] = (
    p_dot_q<p2, p0>(dim_t {}) - p_dot_q<p0, p0>(dim_t {}) +
    dPt_dP<p0, p1, p2, 1>(dim_t {})*lam[0] +
    dPt_dP<p0, p1, p2, 2>(dim_t {})*lam[1]
    // w.dPt_dP[1]*lam[0] +
    // w.dPt_dP[2]*lam[1]
  )/w.l_lam;
  CHECK(w.dPt_nu_lam[1]);
}

template <cost_func F, int n>
void set_lambda(F_wkspc<F, 2> & w, double const * p0, double const * p1,
                double const * p2, double const * lam)
{
  set_lambda_common(w, lam);

  // TODO: it isn't very efficient to recompute these over and over
  // again, when we could just do it once and store the results
  double p0_dot_p0 = dot<n>(p0, p0),
    dp1_dot_p0 = dot<n>(p1, p0) - p0_dot_p0,
    dp2_dot_p0 = dot<n>(p2, p0) - p0_dot_p0;

  w.l_lam = std::sqrt(
    p0_dot_p0 +
    2*(dp1_dot_p0*lam[0] + dp2_dot_p0*lam[1]) +
    w.dPt_dP[0]*lam[0]*lam[0] +
    2*w.dPt_dP[1]*lam[0]*lam[1] +
    w.dPt_dP[2]*lam[1]*lam[1]);
  CHECK(w.l_lam);

  w.dPt_nu_lam[0] = (dp1_dot_p0 + w.dPt_dP[0]*lam[0] + w.dPt_dP[1]*lam[1])/w.l_lam;
  CHECK(w.dPt_nu_lam[0]);

  w.dPt_nu_lam[1] = (dp2_dot_p0 + w.dPt_dP[1]*lam[0] + w.dPt_dP[2]*lam[1])/w.l_lam;
  CHECK(w.dPt_nu_lam[1]);
}

template <cost_func F, int n>
void set_lambda(F0_wkspc<2, fac_wkspc<2>> & w,
                double const * p0, double const * p1, double const * p2,
                double const * p_fac,
                double const * lam)
{
  w.tau_lam = w.tau0 + w.dtau[0]*lam[0] + w.dtau[1]*lam[1];

  // TODO: it isn't very efficient to recompute these over and over
  // again, when we could just do it once and store the results
  double p0_dot_p0 = dot<n>(p0, p0),
    dp1_dot_p0 = dot<n>(p1, p0) - p0_dot_p0,
    dp2_dot_p0 = dot<n>(p2, p0) - p0_dot_p0;

  w.l_lam = std::sqrt(
    p0_dot_p0 +
    2*(dp1_dot_p0*lam[0] + dp2_dot_p0*lam[1]) +
    w.dPt_dP[0]*lam[0]*lam[0] +
    2*w.dPt_dP[1]*lam[0]*lam[1] +
    w.dPt_dP[2]*lam[1]*lam[1]);
  CHECK(w.l_lam);

  w.dPt_nu_lam[0] = (dp1_dot_p0 + w.dPt_dP[0]*lam[0] + w.dPt_dP[1]*lam[1])/w.l_lam;
  CHECK(w.dPt_nu_lam[0]);

  w.dPt_nu_lam[1] = (dp2_dot_p0 + w.dPt_dP[1]*lam[0] + w.dPt_dP[2]*lam[1])/w.l_lam;
  CHECK(w.dPt_nu_lam[1]);

  double p0_dot_p_fac = dot<3>(p0, p_fac),
    dp1_dot_p_fac = dot<n>(p1, p_fac) - p0_dot_p_fac,
    dp2_dot_p_fac = dot<n>(p2, p_fac) - p0_dot_p_fac;

  w.l_fac_lam = w.l_lam*w.l_lam -
    2*(p0_dot_p_fac + lam[0]*dp1_dot_p_fac + lam[1]*dp2_dot_p_fac) +
    dot<n>(p_fac, p_fac);
  if (w.l_fac_lam > 1e1*EPS(double)) {
    w.l_fac_lam = std::sqrt(w.l_fac_lam);
  } else {
    w.l_fac_lam = 0;
  }

  CHECK(w.l_fac_lam);

  if (w.l_fac_lam > 1e1*EPS(double)) {
    w.dPt_nu_fac_lam[0] = (w.l_lam*w.dPt_nu_lam[0] - dp1_dot_p_fac)/w.l_fac_lam;
    CHECK(w.dPt_nu_fac_lam[0]);
    w.dPt_nu_fac_lam[1] = (w.l_lam*w.dPt_nu_lam[1] - dp2_dot_p_fac)/w.l_fac_lam;
    CHECK(w.dPt_nu_fac_lam[1]);
  } else {
    w.dPt_nu_fac_lam[0] = 0;
    w.dPt_nu_fac_lam[1] = 0;
  }
}

template <cost_func F, int n>
void set_lambda(F1_wkspc<2, fac_wkspc<2>> & w,
                double const * p0, double const * p1, double const * p2,
                double const * p_fac,
                double const * lam)
{
  set_lambda<F, n>(static_cast<F0_wkspc<2, fac_wkspc<2>> &>(w), p0, p1, p2, p_fac, lam);

  w.sh_lam = w.sh_bar + w.theta_h_ds[0]*lam[0] + w.theta_h_ds[1]*lam[1];
  CHECK(w.sh_lam);
}

template <int d>
void eval(eval_wkspc<d> const & w, double & f)
{
  f = w.u_lam + w.sh_lam*w.l_lam;
  CHECK(f);
}

template <int d>
void eval(fac_wkspc<d> const & w, double & f)
{
  f = w.tau_lam + w.sh_fac*w.l_fac_lam + w.sh_lam*w.l_lam;
  CHECK(f);
}

template <int d>
void eval_mp1_fix(
  eval_wkspc<d> const & w,
  double s, double s0, double s1, double s2, double h,
  double const * lam, double & f)
{
  f = w.u_lam + h*(s + s0 + (s1 - s0)*lam[0] + (s2 - s0)*lam[1])*w.l_lam/2;
  CHECK(f);
}

template <int d>
void eval_mp1_fix(
  fac_wkspc<d> const & w,
  double s, double s0, double s1, double s2, double h,
  double const * lam, double & f)
{
  f = w.tau_lam + w.sh_fac*w.l_fac_lam +
    h*(s + s0 + (s1 - s0)*lam[0] + (s2 - s0)*lam[1])*w.l_lam/2;
  CHECK(f);
}

template <int d>
void grad(F0_wkspc<d> const & w, double * df)
{
  for (int i = 0; i < d; ++i) df[i] = w.du[i] + w.sh_lam*w.dPt_nu_lam[i];
  for (int i = 0; i < d; ++i) CHECK(df[i]);
}

template <int d>
void grad(F0_fac_wkspc<d> const & w, double * df)
{
  for (int i = 0; i < d; ++i) df[i] = w.dtau[i] + w.sh_lam*w.dPt_nu_lam[i];
  if (w.l_fac_lam > 1e1*EPS(double)) {
    for (int i = 0; i < d; ++i) df[i] += w.sh_fac*w.dPt_nu_fac_lam[i];
  }
  for (int i = 0; i < d; ++i) CHECK(df[i]);
}

template <int d>
void grad(F1_wkspc<d> const & w, double * df)
{
  for (int i = 0; i < d; ++i) {
    df[i] = w.du[i] + w.theta_h_ds[i]*w.l_lam + w.sh_lam*w.dPt_nu_lam[i];
  }
  for (int i = 0; i < d; ++i) CHECK(df[i]);
}

template <int d>
void grad(F1_fac_wkspc<d> const & w, double * df)
{
  for (int i = 0; i < d; ++i) {
    df[i] = w.dtau[i] + w.theta_h_ds[i]*w.l_lam + w.sh_lam*w.dPt_nu_lam[i];
  }
  if (w.l_fac_lam > 1e1*EPS(double)) {
    for (int i = 0; i < d; ++i) df[i] += w.sh_fac*w.dPt_nu_fac_lam[i];
  }
  for (int i = 0; i < d; ++i) CHECK(df[i]);
}

template <int n, int p0, int p1, int p2>
void hess(F0_wkspc<2> const & w, double * d2f)
{
  using namespace bitops;
  using dim_t = dim<n>;

  double const tmp = w.sh_lam/w.l_lam;

  d2f[0] = tmp*(dPt_dP<p0,p1,p2,0>(dim_t {}) - w.dPt_nu_lam[0]*w.dPt_nu_lam[0]);
  CHECK(d2f[0]);

  d2f[1] = tmp*(dPt_dP<p0,p1,p2,1>(dim_t {}) - w.dPt_nu_lam[0]*w.dPt_nu_lam[1]);
  CHECK(d2f[1]);

  d2f[2] = tmp*(dPt_dP<p0,p1,p2,2>(dim_t {}) - w.dPt_nu_lam[1]*w.dPt_nu_lam[1]);
  CHECK(d2f[2]);
}

template <int d>
void hess(F0_wkspc<d> const & w, double * d2f)
{
  double tmp = w.sh_lam/w.l_lam;
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      int k = (d - 1)*i + j;
      d2f[k] = tmp*(w.dPt_dP[k] - w.dPt_nu_lam[i]*w.dPt_nu_lam[j]);
    }
  }
  for (int k = 0; k < sym_mat_size(d); ++k) CHECK(d2f[k]);
}

// TODO: this and the F1 fac hess below are identical, should try to
// combine them into one function
template <int d>
void hess(F0_fac_wkspc<d> const & w, double * d2f)
{
  {
    double tmp = w.sh_lam/w.l_lam;
    for (int i = 0; i < d; ++i) {
      for (int j = i; j < d; ++j) {
        int k = (d - 1)*i + j;
        d2f[k] = tmp*(w.dPt_dP[k] - w.dPt_nu_lam[i]*w.dPt_nu_lam[j]);
      }
    }
  }
  if (w.l_fac_lam > EPS(double)) {
    double tmp = w.sh_fac/w.l_fac_lam;
    for (int i = 0; i < d; ++i) {
      for (int j = i; j < d; ++j) {
        int k = (d - 1)*i + j;
        d2f[k] += tmp*(w.dPt_dP[k] - w.dPt_nu_fac_lam[i]*w.dPt_nu_fac_lam[j]);
      }
    }
  }
  for (int k = 0; k < sym_mat_size(d); ++k) CHECK(d2f[k]);
}

template <int n, int p0, int p1, int p2>
void hess(F1_wkspc<2> const & w, double * d2f)
{
  hess<n, p0, p1, p2>(static_cast<F0_wkspc<2> const &>(w), d2f);

  d2f[0] += 2*w.theta_h_ds[0]*w.dPt_nu_lam[0];
  d2f[1] += w.theta_h_ds[0]*w.dPt_nu_lam[1] + w.theta_h_ds[1]*w.dPt_nu_lam[0];
  d2f[2] += 2*w.theta_h_ds[1]*w.dPt_nu_lam[1];
}

template <int d>
void hess(F1_wkspc<d> const & w, double * d2f)
{
  double tmp = w.sh_lam/w.l_lam;
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      int k = (d - 1)*i + j;
      d2f[k] = w.theta_h_ds[i]*w.dPt_nu_lam[j] +
        w.theta_h_ds[j]*w.dPt_nu_lam[i] +
        tmp*(w.dPt_dP[k] - w.dPt_nu_lam[i]*w.dPt_nu_lam[j]);
    }
  }
  for (int k = 0; k < sym_mat_size(d); ++k) CHECK(d2f[k]);
}

template <int d>
void hess(F1_fac_wkspc<d> const & w, double * d2f)
{
  {
    double tmp = w.sh_lam/w.l_lam;
    for (int i = 0; i < d; ++i) {
      for (int j = i; j < d; ++j) {
        int k = (d - 1)*i + j;
        d2f[k] = w.theta_h_ds[i]*w.dPt_nu_lam[j] +
          w.theta_h_ds[j]*w.dPt_nu_lam[i] +
          tmp*(w.dPt_dP[k] - w.dPt_nu_lam[i]*w.dPt_nu_lam[j]);
      }
    }
  }
  if (w.l_fac_lam > EPS(double)) {
    double tmp = w.sh_fac/w.l_fac_lam;
    for (int i = 0; i < d; ++i) {
      for (int j = i; j < d; ++j) {
        int k = (d - 1)*i + j;
        d2f[k] += tmp*(w.dPt_dP[k] - w.dPt_nu_fac_lam[i]*w.dPt_nu_fac_lam[j]);
      }
    }
  }
  for (int k = 0; k < sym_mat_size(d); ++k) CHECK(d2f[k]);
}

// TODO: it would be nice to get rid of these and be able to just pass
// these functions to SQP directly (for example, as template parameters)

template <cost_func F, int... ps>
struct cost_functor_bv;

template <cost_func F, int n, int p0, int p1, int p2>
struct cost_functor_bv<F, n, p0, p1, p2>
{
  using wkspc = F_wkspc<F, 2>;
  cost_functor_bv(F_wkspc<F, 2> & w): w {w} {}
  inline void set_lambda(double const * lam) {
    ::set_lambda<F, n, p0, p1, p2>(w, lam);
  }
  inline void eval(double & f) const {::eval(w, f);}
  inline void grad(double * df) const {::grad(w, df);}
  inline void hess(double * d2f) const {::hess<n, p0, p1, p2>(w, d2f);}
  wkspc & w;
};

template <cost_func F, int n>
struct cost_functor
{
  using wkspc = F_wkspc<F, 2>;
  cost_functor(F_wkspc<F, 2> & w,
               double const * p0, double const * p1, double const * p2):
    w {w}, p0 {p0}, p1 {p1}, p2 {p2} {}
  inline void set_lambda(double const * lam) {
    ::set_lambda<F, n>(w, p0, p1, p2, lam);
  }
  inline void eval(double & f) const {::eval(w, f);}
  inline void grad(double * df) const  {::grad(w, df);}
  inline void hess(double * d2f) const {::hess(w, d2f);}
  wkspc & w;
  double const * p0, * p1, * p2;
};

template <cost_func F, int n>
struct cost_functor_fac
{
  using wkspc = F_fac_wkspc<F, 2>;
  cost_functor_fac(F_fac_wkspc<F, 2> & w,
                   double const * p0, double const * p1, double const * p2,
                   double const * p_fac):
    w {w}, p0 {p0}, p1 {p1}, p2 {p2}, p_fac {p_fac} {}
  inline void set_lambda(double const * lam) {
    ::set_lambda<F, n>(w, p0, p1, p2, p_fac, lam);
  }
  inline void eval(double & f) const {::eval(w, f);}
  inline void grad(double * df) const {::grad(w, df);}
  inline void hess(double * d2f) const {::hess(w, d2f);}
  wkspc & w;
  double const * p0, * p1, * p2, * p_fac;
};


#endif // __COST_FUNCS_HPP__
