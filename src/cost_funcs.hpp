#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

#include <assert.h>
#include <math.h>

// TODO: remove me
#include <type_traits>

#include "bitops.hpp"
#include "common.hpp"
#include "vecmath.hpp"

enum class cost_func {mp0, mp1, rhr};

constexpr auto MP0 = cost_func::mp0;
constexpr auto MP1 = cost_func::mp1;
constexpr auto RHR = cost_func::rhr;

constexpr int sym_mat_size(int d) {
  return ((d + 1)*d)/2;
}

template <int d> inline void check_lambda(double const * lam);

template <>
inline void check_lambda<2>(double const * lam) {
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(lam[0] >= -eps<double>);
  assert(lam[1] >= -eps<double>);
  assert(lam[0] + lam[1] <= 1 + eps<double>);
#else
  (void) lam;
#endif
}

inline void check(double x) {
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(!isinf(x));
  assert(!isnan(x));
#else
  (void) x;
#endif
}

template <int d>
void lagmults(double const * lam, double const * df, double const * d2f,
              double * mu, int * k);

template <int n, int d> struct qr_wkspc {};

template <> struct qr_wkspc<3, 2> {
  double q1[3];
  double q2[3];
  double r[3];
  double Qt_p0[2];
  double numer;

  void init(double const * p0, double const * p1, double const * p2) {
    // Compute 3x2 reduced QR decomposition of [p1 - p0, p2 - p0].
    sub<3>(p1, p0, q1);
    r[0] = norm2<3>(q1);
    scal_inplace<3>(1/r[0], q1);
    sub<3>(p2, p0, q2);
    r[1] = dot<3>(q2, q1);
    axpy<3>(-r[1], q1, q2, q2);
    r[2] = norm2<3>(q2);
    scal_inplace<3>(1/r[2], q2);

    // Compute entries of Q'*p0:
    Qt_p0[0] = dot<3>(q1, p0);
    Qt_p0[1] = dot<3>(q2, p0);

    // Compute the numerator in the square root for the exact
    // solution (TODO: add reference to equation in paper).
    numer = dot<3>(p0, p0) - dot<2>(Qt_p0, Qt_p0);
  }
};

template <int d> struct geom_wkspc {};

template <>
struct geom_wkspc<2> {
  double p0t_p0, dPt_p0[2], dPt_dP[sym_mat_size(2)];

  template <int n>
  void init(double const * p0, double const * p1, double const * p2) {
    double dp1[n], dp2[n];
    sub<n>(p1, p0, dp1);
    sub<n>(p2, p0, dp2);

    p0t_p0 = dot<n>(p0, p0);
    dPt_p0[0] = dot<n>(dp1, p0);
    dPt_p0[1] = dot<n>(dp2, p0);
    dPt_dP[0] = dot<n>(dp1, dp1);
    dPt_dP[1] = dot<n>(dp1, dp2);
    dPt_dP[2] = dot<n>(dp2, dp2);
  }
};

template <int d> struct geom_fac_wkspc {};

template <>
struct geom_fac_wkspc<2>: geom_wkspc<2> {
  double pft_pf, dPt_pf[2], pft_p0, L0, dL[2];

  template <int n>
  void init(double const * p0, double const * p1, double const * p2,
            double const * pf) {
    double dp1[n], dp2[n];
    sub<n>(p1, p0, dp1);
    sub<n>(p2, p0, dp2);

    this->p0t_p0 = dot<n>(p0, p0);
    this->dPt_p0[0] = dot<n>(dp1, p0);
    this->dPt_p0[1] = dot<n>(dp2, p0);
    this->dPt_dP[0] = dot<n>(dp1, dp1);
    this->dPt_dP[1] = dot<n>(dp1, dp2);
    this->dPt_dP[2] = dot<n>(dp2, dp2);

    pft_pf = dot<n>(pf, pf);
    dPt_pf[0] = dot<n>(dp1, pf);
    dPt_pf[1] = dot<n>(dp2, pf);
    pft_p0 = dot<n>(pf, p0);
    L0 = dist2<n>(pf, p0);
    dL[0] = dist2<n>(pf, p1) - L0;
    dL[1] = dist2<n>(pf, p2) - L0;
  }
};

template <int d>
struct base_wkspc {
  double sh_lam, l_lam, dPt_nu_lam[d];
};

template <int d>
struct F0_wkspc: base_wkspc<d> {
  double u0, du[d], u_lam;
};

template <int d>
struct fac_wkspc: base_wkspc<d> {
  double tau0, dtau[d], tau_lam, sh_fac, l_fac_lam, dPt_nu_fac_lam[d];
};

template <int d, class base = F0_wkspc<d>>
struct F1_wkspc: base {
  double sh_bar, theta_h_ds[d];
};

template <cost_func F, int d>
using F_wkspc = std::conditional_t<F == MP1, F1_wkspc<d>, F0_wkspc<d>>;

template <int d>
using F0_fac_wkspc = fac_wkspc<d>;

template <int d>
using F1_fac_wkspc = F1_wkspc<d, fac_wkspc<d>>;

template <cost_func F, int d>
using F_fac_wkspc = std::conditional_t<
  F == MP1, F1_fac_wkspc<d>, F0_fac_wkspc<d>>;

inline void check_args(double u0, double u1, double u2, double s,
                       double s0, double s1, double s2, double h)
{
#if EIKONAL_DEBUG && !RELWITHDEBINFO
  assert(u0 >= 0);
  assert(u1 >= 0);
  assert(u2 >= 0);
  assert(s >= 0);
  assert(s0 >= 0);
  assert(s1 >= 0);
  assert(s2 >= 0);
  assert(h > 0);
#else
  (void) u0;
  (void) u1;
  (void) u2;
  (void) s;
  (void) s0;
  (void) s1;
  (void) s2;
  (void) h;
#endif
}

template <cost_func F>
void set_args(F0_wkspc<2> & w,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h)
{
  check_args(u0, u1, u2, s, s0, s1, s2, h);

  w.u0 = u0;
  w.du[0] = u1 - u0;
  w.du[1] = u2 - u0;

  w.sh_lam = (F == RHR ? s : (s + (s0 + s1 + s2)/3)/2)*h;
}

template <cost_func F>
void set_args(F1_wkspc<2> & w,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h)
{
  static_assert(F == MP1, "Cost function must be MP1");

  check_args(u0, u1, u2, s, s0, s1, s2, h);

  set_args<F>(static_cast<F0_wkspc<2> &>(w), u0, u1, u2, s, s0, s1, s2, h);

  w.sh_bar = (s + s0)*h/2;

  w.theta_h_ds[0] = h*(s1 - s0)/2;
  w.theta_h_ds[1] = h*(s2 - s0)/2;
}

template <cost_func F>
void set_args(F0_fac_wkspc<2> & w, geom_fac_wkspc<2> const & g,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h, double s_fac)
{
  check_args(u0, u1, u2, s, s0, s1, s2, h);
  assert(s_fac >= 0);

  w.sh_fac = s_fac*h;

  w.tau0 = u0 - w.sh_fac*g.L0;
  w.dtau[0] = u1 - u0 - w.sh_fac*g.dL[0];
  w.dtau[1] = u2 - u0 - w.sh_fac*g.dL[1];

  w.sh_lam = (F == RHR ? s : (s + (s0 + s1 + s2)/3)/2)*h;
}

template <cost_func F>
void set_args(F1_fac_wkspc<2> & w, geom_fac_wkspc<2> const & g,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h, double s_fac)
{
  static_assert(F == MP1, "Cost function must be MP1");

  check_args(u0, u1, u2, s, s0, s1, s2, h);
  assert(s_fac >= 0);

  set_args<F>(
    static_cast<F0_fac_wkspc<2> &>(w), g,
    u0, u1, u2, s, s0, s1, s2, h, s_fac);

  w.sh_bar = (s + s0)*h/2;

  w.theta_h_ds[0] = h*(s1 - s0)/2;
  w.theta_h_ds[1] = h*(s2 - s0)/2;
}

inline void set_lambda_common(F0_wkspc<2> & w, double const * lam)
{
  check_lambda<2>(lam);

  w.u_lam = w.u0 + w.du[0]*lam[0] + w.du[1]*lam[1];
}

inline void set_lambda_common(F1_wkspc<2> & w, double const * lam)
{
  check_lambda<2>(lam);

  set_lambda_common(static_cast<F0_wkspc<2> &>(w), lam);

  w.sh_lam = w.sh_bar + w.theta_h_ds[0]*lam[0] + w.theta_h_ds[1]*lam[1];
}

#define __p0t_p0 p_dot_q<p0, p0>(dim_t {})
#define __dPt_p0(i) p_dot_q<p##i, p0>(dim_t {}) - p_dot_q<p0, p0>(dim_t {})
#define __dPt_dP(i) dPt_dP<p0, p1, p2, i>(dim_t {})

template <cost_func F, int n, int p0, int p1, int p2>
void set_lambda(F_wkspc<F, 2> & w, double const * lam)
{
  using namespace bitops;
  using dim_t = dim<n>;

  check_lambda<2>(lam);

  set_lambda_common(w, lam);

  w.dPt_nu_lam[0] = __dPt_dP(0)*lam[0] + __dPt_dP(1)*lam[1] + __dPt_p0(1);
  w.dPt_nu_lam[1] = __dPt_dP(1)*lam[0] + __dPt_dP(2)*lam[1] + __dPt_p0(2);

  w.l_lam = sqrt(
    __p0t_p0 +
    (__dPt_p0(1) + w.dPt_nu_lam[0])*lam[0] +
    (__dPt_p0(2) + w.dPt_nu_lam[1])*lam[1]
    );
  check(w.l_lam);

  w.dPt_nu_lam[0] /= w.l_lam;
  w.dPt_nu_lam[1] /= w.l_lam;

  check(w.dPt_nu_lam[0]);
  check(w.dPt_nu_lam[1]);
}

#undef __p0t_p0
#undef __dPt_p0
#undef __dPt_dP

template <cost_func F>
void set_lambda(F_wkspc<F, 2> & w, geom_wkspc<2> & g, double const * lam)
{
  check_lambda<2>(lam);

  set_lambda_common(w, lam);

  w.dPt_nu_lam[0] = g.dPt_dP[0]*lam[0] + g.dPt_dP[1]*lam[1] + g.dPt_p0[0];
  w.dPt_nu_lam[1] = g.dPt_dP[1]*lam[0] + g.dPt_dP[2]*lam[1] + g.dPt_p0[1];

  w.l_lam = sqrt(
    g.p0t_p0 +
    (g.dPt_p0[0] + w.dPt_nu_lam[0])*lam[0] +
    (g.dPt_p0[1] + w.dPt_nu_lam[1])*lam[1]
  );
  check(w.l_lam);

  w.dPt_nu_lam[0] /= w.l_lam;
  w.dPt_nu_lam[1] /= w.l_lam;

  check(w.dPt_nu_lam[0]);
  check(w.dPt_nu_lam[1]);
}

template <cost_func F>
void set_lambda(F0_fac_wkspc<2> & w, geom_fac_wkspc<2> const & g,
                double const * lam)
{
  check_lambda<2>(lam);

  w.tau_lam = w.tau0 + w.dtau[0]*lam[0] + w.dtau[1]*lam[1];

  // TODO: we can probably simply this quite a lot, following the same
  // idea in the unfactored `set_lambda' function above.

  w.l_lam = sqrt(
    g.p0t_p0 +
    2*(g.dPt_p0[0]*lam[0] + g.dPt_p0[1]*lam[1]) +
    g.dPt_dP[0]*lam[0]*lam[0] +
    2*g.dPt_dP[1]*lam[0]*lam[1] +
    g.dPt_dP[2]*lam[1]*lam[1]);
  check(w.l_lam);

  w.dPt_nu_lam[0] = (g.dPt_p0[0] + g.dPt_dP[0]*lam[0] + g.dPt_dP[1]*lam[1])/w.l_lam;
  check(w.dPt_nu_lam[0]);

  w.dPt_nu_lam[1] = (g.dPt_p0[1] + g.dPt_dP[1]*lam[0] + g.dPt_dP[2]*lam[1])/w.l_lam;
  check(w.dPt_nu_lam[1]);

  // TODO: the line below with the funny-looking check before we take
  // the square root is a hack to get around the fact that if we take
  // the square root of a very small number, like 1e-16, then the
  // result is 1e-8. This can cause problems when p_fac corresponds
  // with one of the update simplices (so, near the point source). The
  // result can be needlessly inflated numerical error that propagates
  // throughout the entire solution.
  w.l_fac_lam = w.l_lam*w.l_lam - // TODO: set l_fac_lam to
                                  // unnormalize l_lam first before
                                  // adding the rest to it and
                                  // normalizing l_lam
    2*(g.pft_p0 + lam[0]*g.dPt_pf[0] + lam[1]*g.dPt_pf[1]) + g.pft_pf;
  w.l_fac_lam = w.l_fac_lam > 1e1*eps<double> ? sqrt(w.l_fac_lam) : 0;
  check(w.l_fac_lam);

  w.dPt_nu_fac_lam[0] = (w.l_lam*w.dPt_nu_lam[0] - g.dPt_pf[0])/w.l_fac_lam;
  w.dPt_nu_fac_lam[1] = (w.l_lam*w.dPt_nu_lam[1] - g.dPt_pf[1])/w.l_fac_lam;
  if (!isinf(w.dPt_nu_fac_lam[0])) {
    w.dPt_nu_fac_lam[0] = w.dPt_nu_fac_lam[1] = 0;
  }
}

template <cost_func F>
void set_lambda(F1_fac_wkspc<2> & w, geom_fac_wkspc<2> const & g,
                double const * lam)
{
  check_lambda<2>(lam);

  set_lambda<F>(static_cast<F0_fac_wkspc<2> &>(w), g, lam);

  w.sh_lam = w.sh_bar + w.theta_h_ds[0]*lam[0] + w.theta_h_ds[1]*lam[1];
  check(w.sh_lam);
}

template <int d>
void eval(F0_wkspc<d> const & w, double & f)
{
  f = w.u_lam + w.sh_lam*w.l_lam;
  check(f);
}

template <int d>
void eval(fac_wkspc<d> const & w, double & f)
{
  f = w.tau_lam + w.sh_fac*w.l_fac_lam + w.sh_lam*w.l_lam;
  check(f);
}

inline void grad(F0_wkspc<2> const & w, double * df)
{
  df[0] = w.du[0] + w.sh_lam*w.dPt_nu_lam[0];
  df[1] = w.du[1] + w.sh_lam*w.dPt_nu_lam[1];

  check(df[0]);
  check(df[1]);
}

inline void grad(F0_fac_wkspc<2> const & w, double * df)
{
  df[0] = w.dtau[0] + w.sh_lam*w.dPt_nu_lam[0];
  df[1] = w.dtau[1] + w.sh_lam*w.dPt_nu_lam[1];

  if (w.l_fac_lam > 1e1*eps<double>) {
    df[0] += w.sh_fac*w.dPt_nu_fac_lam[0];
    df[1] += w.sh_fac*w.dPt_nu_fac_lam[1];
  }

  check(df[0]);
  check(df[1]);
}

inline void grad(F1_wkspc<2> const & w, double * df)
{
  df[0] = w.du[0] + w.theta_h_ds[0]*w.l_lam + w.sh_lam*w.dPt_nu_lam[0];
  df[1] = w.du[1] + w.theta_h_ds[1]*w.l_lam + w.sh_lam*w.dPt_nu_lam[1];

  check(df[0]);
  check(df[1]);
}

inline void grad(F1_fac_wkspc<2> const & w, double * df)
{
  df[0] = w.dtau[0] + w.theta_h_ds[0]*w.l_lam + w.sh_lam*w.dPt_nu_lam[0];
  df[1] = w.dtau[1] + w.theta_h_ds[1]*w.l_lam + w.sh_lam*w.dPt_nu_lam[1];

  if (w.l_fac_lam > 1e1*eps<double>) {
    df[0] += w.sh_fac*w.dPt_nu_fac_lam[0];
    df[1] += w.sh_fac*w.dPt_nu_fac_lam[1];
  }

  check(df[0]);
  check(df[1]);
}

#define __dPt_dP(i) dPt_dP<p0, p1, p2, i>(dim_t {})

template <int n, int p0, int p1, int p2>
void hess(F0_wkspc<2> const & w, double * d2f)
{
  using namespace bitops;
  using dim_t = dim<n>;

  double const tmp = w.sh_lam/w.l_lam;

  d2f[0] = tmp*(__dPt_dP(0) - w.dPt_nu_lam[0]*w.dPt_nu_lam[0]);
  d2f[1] = tmp*(__dPt_dP(1) - w.dPt_nu_lam[0]*w.dPt_nu_lam[1]);
  d2f[2] = tmp*(__dPt_dP(2) - w.dPt_nu_lam[1]*w.dPt_nu_lam[1]);

  check(d2f[0]);
  check(d2f[1]);
  check(d2f[2]);
}

#undef __dPt_dP

inline void hess(F0_wkspc<2> const & w, geom_wkspc<2> const & g, double * d2f)
{
  double tmp = w.sh_lam/w.l_lam;

  d2f[0] = tmp*(g.dPt_dP[0] - w.dPt_nu_lam[0]*w.dPt_nu_lam[0]);
  d2f[1] = tmp*(g.dPt_dP[1] - w.dPt_nu_lam[0]*w.dPt_nu_lam[1]);
  d2f[2] = tmp*(g.dPt_dP[2] - w.dPt_nu_lam[1]*w.dPt_nu_lam[1]);

  check(d2f[0]);
  check(d2f[1]);
  check(d2f[2]);
}

inline void hess(F0_fac_wkspc<2> const & w,
                 geom_fac_wkspc<2> const & g,
                 double * d2f)
{
  double tmp = w.sh_lam/w.l_lam;
  d2f[0] = tmp*(g.dPt_dP[0] - w.dPt_nu_lam[0]*w.dPt_nu_lam[0]);
  d2f[1] = tmp*(g.dPt_dP[1] - w.dPt_nu_lam[0]*w.dPt_nu_lam[1]);
  d2f[2] = tmp*(g.dPt_dP[2] - w.dPt_nu_lam[1]*w.dPt_nu_lam[1]);

  tmp = w.sh_fac/w.l_fac_lam;
  if (!isinf(tmp)) {
    d2f[0] += tmp*(g.dPt_dP[0] - w.dPt_nu_fac_lam[0]*w.dPt_nu_fac_lam[0]);
    d2f[1] += tmp*(g.dPt_dP[1] - w.dPt_nu_fac_lam[0]*w.dPt_nu_fac_lam[1]);
    d2f[2] += tmp*(g.dPt_dP[2] - w.dPt_nu_fac_lam[1]*w.dPt_nu_fac_lam[1]);
  }

  check(d2f[0]);
  check(d2f[1]);
  check(d2f[2]);
}

template <int n, int p0, int p1, int p2>
void hess(F1_wkspc<2> const & w, double * d2f)
{
  hess<n, p0, p1, p2>(static_cast<F0_wkspc<2> const &>(w), d2f);

  d2f[0] += 2*w.theta_h_ds[0]*w.dPt_nu_lam[0];
  d2f[1] += w.theta_h_ds[0]*w.dPt_nu_lam[1] + w.theta_h_ds[1]*w.dPt_nu_lam[0];
  d2f[2] += 2*w.theta_h_ds[1]*w.dPt_nu_lam[1];
}

inline void hess(F1_wkspc<2> const & w, geom_wkspc<2> const & g, double * d2f)
{
  double tmp = w.sh_lam/w.l_lam;

  d2f[0] = 2*w.theta_h_ds[0]*w.dPt_nu_lam[0] +
    tmp*(g.dPt_dP[0] - w.dPt_nu_lam[0]*w.dPt_nu_lam[0]);

  d2f[1] = w.theta_h_ds[0]*w.dPt_nu_lam[1] + w.theta_h_ds[1]*w.dPt_nu_lam[0] +
    tmp*(g.dPt_dP[1] - w.dPt_nu_lam[0]*w.dPt_nu_lam[1]);

  d2f[2] = 2*w.theta_h_ds[1]*w.dPt_nu_lam[1] +
    tmp*(g.dPt_dP[2] - w.dPt_nu_lam[1]*w.dPt_nu_lam[1]);

  check(d2f[0]);
  check(d2f[1]);
  check(d2f[2]);
}

inline void hess(F1_fac_wkspc<2> const & w,
                 geom_fac_wkspc<2> const & g,
                 double * d2f)
{
  double tmp = w.sh_lam/w.l_lam;

  d2f[0] = 2*w.theta_h_ds[0]*w.dPt_nu_lam[0] +
    tmp*(g.dPt_dP[0] - w.dPt_nu_lam[0]*w.dPt_nu_lam[0]);

  d2f[1] = w.theta_h_ds[0]*w.dPt_nu_lam[1] + w.theta_h_ds[1]*w.dPt_nu_lam[0] +
    tmp*(g.dPt_dP[1] - w.dPt_nu_lam[0]*w.dPt_nu_lam[1]);

  d2f[2] = 2*w.theta_h_ds[1]*w.dPt_nu_lam[1] +
    tmp*(g.dPt_dP[2] - w.dPt_nu_lam[1]*w.dPt_nu_lam[1]);

  if (w.l_fac_lam > 1e1*eps<double>) {
    tmp = w.sh_fac/w.l_fac_lam;
    d2f[0] += tmp*(g.dPt_dP[0] - w.dPt_nu_fac_lam[0]*w.dPt_nu_fac_lam[0]);
    d2f[1] += tmp*(g.dPt_dP[1] - w.dPt_nu_fac_lam[0]*w.dPt_nu_fac_lam[1]);
    d2f[2] += tmp*(g.dPt_dP[2] - w.dPt_nu_fac_lam[1]*w.dPt_nu_fac_lam[1]);
  }

  check(d2f[0]);
  check(d2f[1]);
  check(d2f[2]);
}

// TODO: it would be nice to get rid of these and be able to just pass
// these functions to SQP directly (for example, as template parameters)

template <cost_func F, int... ps>
struct cost_functor_bv;

template <cost_func F, int n, int p0, int p1, int p2>
struct cost_functor_bv<F, n, p0, p1, p2>
{
  cost_functor_bv(F_wkspc<F, 2> & w): w {w} {}
  inline void set_lambda(double const * lam) { ::set_lambda<F, n, p0, p1, p2>(w, lam); }
  inline void eval(double & f) const {::eval(w, f);}
  inline void grad(double * df) const {::grad(w, df);}
  inline void hess(double * d2f) const {::hess<n, p0, p1, p2>(w, d2f);}
  F_wkspc<F, 2> & w;
};

template <cost_func F, int n, int d>
struct cost_functor
{
  cost_functor(F_wkspc<F, d> & w, geom_wkspc<d> & g): w {w}, g {g} {}
  inline void set_lambda(double const * lam) {::set_lambda<F>(w, g, lam);}
  inline void eval(double & f) const {::eval(w, f);}
  inline void grad(double * df) const  {::grad(w, df);}
  inline void hess(double * d2f) const {::hess(w, g, d2f);}
  F_wkspc<F, d> & w;
  geom_wkspc<d> & g;
  qr_wkspc<n, d> * qr {nullptr};
};

template <cost_func F, int n, int d>
struct cost_functor_fac
{
  cost_functor_fac(F_fac_wkspc<F, d> & w, geom_fac_wkspc<d> & g): w {w}, g {g} {}
  inline void set_lambda(double const * lam) {::set_lambda<F>(w, g, lam); }
  inline void eval(double & f) const {::eval(w, f);}
  inline void grad(double * df) const {::grad(w, df);}
  inline void hess(double * d2f) const {::hess(w, g, d2f);}
  F_fac_wkspc<F, d> & w;
  geom_fac_wkspc<d> & g;
};

inline void eval_mp1_fix(F0_wkspc<2> const & w,
                  double s, double s0, double s1, double s2, double h,
                  double const * lam, double & f)
{
  f = w.u_lam + h*(s + s0 + (s1 - s0)*lam[0] + (s2 - s0)*lam[1])*w.l_lam/2;
  check(f);
}

inline void eval_mp1_fix(fac_wkspc<2> const & w,
                  double s, double s0, double s1, double s2, double h,
                  double const * lam, double & f)
{
  f = w.tau_lam + w.sh_fac*w.l_fac_lam +
    h*(s + s0 + (s1 - s0)*lam[0] + (s2 - s0)*lam[1])*w.l_lam/2;
  check(f);
}

#endif // __COST_FUNCS_HPP__
