#ifndef __COST_FUNCS_IMPL_HPP__
#define __COST_FUNCS_IMPL_HPP__

#include <cmath>

#include "update_rules.utils.hpp"

template <char pj, char p0>
constexpr char dp(char i) {
  return component(pj, i) - component(p0, i);
}

#define __dP(i, j) dp<p##i, p0>(j)

#define __dP_dot_prod_n3(i, x)                          \
  (__dP(i, 0)*x[0] + __dP(i, 1)*x[1] + __dP(i, 2)*x[2])

template <char p0, char p1, char p2>
void
F0_bv<p0, p1, p2, 2>::eval_impl(double & f) const
{
  f = _u_lam + _sh*_l;
}

template <char p0, char p1, char p2>
void
F0_bv<p0, p1, p2, 2>::grad_impl(double df[2]) const
{
  df[0] = _du[0] + _sh*__dP_dot_prod_n3(1, _p_lam)/_l;
  df[1] = _du[1] + _sh*__dP_dot_prod_n3(2, _p_lam)/_l;
}

template <char p0, char p1, char p2>
void
F0_bv<p0, p1, p2, 2>::hess_impl(double d2f[2][2]) const
{
  double tmp[6]; // workspace
  tmp[0] = 1 - _p_lam[0]*_p_lam[0]/_q; // (0, 0)
  tmp[1] = -_p_lam[0]*_p_lam[1]/_q;    // (0, 1) & (1, 0)
  tmp[2] = -_p_lam[0]*_p_lam[2]/_q;    // (0, 2) & (2, 0)
  tmp[3] = 1 - _p_lam[1]*_p_lam[1]/_q; // (1, 1)
  tmp[4] = -_p_lam[1]*_p_lam[2]/_q;    // (1, 2) & (2, 1)
  tmp[5] = 1 - _p_lam[2]*_p_lam[2]/_q; // (2, 2)

  double tmp2[2][3];
  // tmp2[0][0] = tmp[0]*__dP(1, 0) + tmp[1]*__dP(1, 1) + tmp[2]*__dP(1, 2);
  tmp2[0][0] = __dP_dot_prod_n3(1, tmp);
  tmp2[0][1] = tmp[1]*__dP(1, 0) + tmp[3]*__dP(1, 1) + tmp[4]*__dP(1, 2);
  tmp2[0][2] = tmp[2]*__dP(1, 0) + tmp[4]*__dP(1, 1) + tmp[5]*__dP(1, 2);
  // tmp2[1][0] = tmp[0]*__dP(2, 0) + tmp[1]*__dP(2, 1) + tmp[2]*__dP(2, 2);
  tmp2[1][0] = __dP_dot_prod_n3(2, tmp);
  tmp2[1][1] = tmp[1]*__dP(2, 0) + tmp[3]*__dP(2, 1) + tmp[4]*__dP(2, 2);
  tmp2[1][2] = tmp[2]*__dP(2, 0) + tmp[4]*__dP(2, 1) + tmp[5]*__dP(2, 2);

  // tmp[0] = __dP(1, 0)*tmp2[0][0] + __dP(1, 1)*tmp2[0][1] + __dP(1, 2)*tmp2[0][2];
  // tmp[1] = __dP(1, 0)*tmp2[1][0] + __dP(1, 1)*tmp2[1][1] + __dP(1, 2)*tmp2[1][2];
  // tmp[2] = __dP(2, 0]*tmp2[1][0] + __dP(2, 1)*tmp2[1][1] + __dP(2, 2)*tmp2[1][2];
  tmp[0] = __dP_dot_prod_n3(1, tmp2[0]);
  tmp[1] = __dP_dot_prod_n3(1, tmp2[1]);
  tmp[2] = __dP_dot_prod_n3(2, tmp2[1]);

  d2f[0][0] = _sh*tmp[0]/_l;
  d2f[0][1] = d2f[1][0] = _sh*tmp[1]/_l;
  d2f[1][1] = _sh*tmp[2]/_l;
}

template <char p0, char p1, char p2>
void
F0_bv<p0, p1, p2, 2>::set_lambda_impl(double const lam[2])
{
  _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];

  _p_lam[0] = component(p0, 0) + lam[0]*__dP(1, 0) + lam[1]*__dP(2, 0);
  _p_lam[1] = component(p0, 1) + lam[0]*__dP(1, 1) + lam[1]*__dP(2, 1);
  _p_lam[2] = component(p0, 2) + lam[0]*__dP(1, 2) + lam[1]*__dP(2, 2);

  // speed this up by rewriting p^T*p?
  _q = _p_lam[0]*_p_lam[0] + _p_lam[1]*_p_lam[1] + _p_lam[2]*_p_lam[2];

  // store _sh/_l instead?
  _l = sqrt(_q);
}

template <char p0, char p1, char p2>
void
F0_bv<p0, p1, p2, 2>::set_args_impl(double const u[3], double s_hat,
                                    double const s[3])
{
  _sh = _h*((1 - _theta)*s_hat + _theta*(s[0] + s[1] + s[2])/3);
  _u0 = u[0];
  _du[0] = u[1] - _u0;
  _du[1] = u[2] - _u0;
}

template <char p0, char p1, char p2>
void
F1_bv<p0, p1, p2, 2>::eval_impl(double & f) const
{
  f = _u_lam + _h*_stheta*_l;
}

template <char p0, char p1, char p2>
void
F1_bv<p0, p1, p2, 2>::grad_impl(double df[2]) const
{
  df[0] = _du[0] + _h*(_theta*_q*_ds[0] + _stheta*__dP_dot_prod_n3(1, _p_lam))/_l;
  df[1] = _du[1] + _h*(_theta*_q*_ds[1] + _stheta*__dP_dot_prod_n3(2, _p_lam))/_l;
}

template <char p0, char p1, char p2>
void
F1_bv<p0, p1, p2, 2>::hess_impl(double d2f[2][2]) const
{
  double tmp[6]; // workspace
  tmp[0] = 1 - _p_lam[0]*_p_lam[0]/_q; // (0, 0)
  tmp[1] = -_p_lam[0]*_p_lam[1]/_q;    // (0, 1) & (1, 0)
  tmp[2] = -_p_lam[0]*_p_lam[2]/_q;    // (0, 2) & (2, 0)
  tmp[3] = 1 - _p_lam[1]*_p_lam[1]/_q; // (1, 1)
  tmp[4] = -_p_lam[1]*_p_lam[2]/_q;    // (1, 2) & (2, 1)
  tmp[5] = 1 - _p_lam[2]*_p_lam[2]/_q; // (2, 2)

  double tmp2[2][3];
  tmp2[0][0] = tmp[0]*__dP(1, 0) + tmp[1]*__dP(1, 1) + tmp[2]*__dP(1, 2);
  tmp2[0][1] = tmp[1]*__dP(1, 0) + tmp[3]*__dP(1, 1) + tmp[4]*__dP(1, 2);
  tmp2[0][2] = tmp[2]*__dP(1, 0) + tmp[4]*__dP(1, 1) + tmp[5]*__dP(1, 2);
  tmp2[1][0] = tmp[0]*__dP(2, 0) + tmp[1]*__dP(2, 1) + tmp[2]*__dP(2, 2);
  tmp2[1][1] = tmp[1]*__dP(2, 0) + tmp[3]*__dP(2, 1) + tmp[4]*__dP(2, 2);
  tmp2[1][2] = tmp[2]*__dP(2, 0) + tmp[4]*__dP(2, 1) + tmp[5]*__dP(2, 2);

  // compute dP'*(I - p*p'/q)*dP
  tmp[0] = __dP(1, 0)*tmp2[0][0] + __dP(1, 1)*tmp2[0][1] + __dP(1, 2)*tmp2[0][2];
  tmp[1] = __dP(1, 0)*tmp2[1][0] + __dP(1, 1)*tmp2[1][1] + __dP(1, 2)*tmp2[1][2];
  tmp[2] = __dP(2, 0)*tmp2[1][0] + __dP(2, 1)*tmp2[1][1] + __dP(2, 2)*tmp2[1][2];

  // compute dP'*p
  tmp[3] = __dP(1, 0)*_p_lam[0] + __dP(1, 1)*_p_lam[1] + __dP(1, 2)*_p_lam[2];
  tmp[4] = __dP(2, 0)*_p_lam[0] + __dP(2, 1)*_p_lam[1] + __dP(2, 2)*_p_lam[2];

  // compute theta*(dP'*p*ds' + ds*p'*dP)
  tmp2[0][0] = 2*_theta*tmp[3]*_ds[0];
  tmp2[0][1] = _theta*(tmp[3]*_ds[1] + tmp[4]*_ds[0]);
  tmp2[0][2] = 2*_theta*tmp[4]*_ds[1];

  // compute Hessian
  d2f[0][0] = _h*(tmp2[0][0] + _stheta*tmp[0])/_l;
  d2f[0][1] = d2f[1][0] = _h*(tmp2[0][1] + _stheta*tmp[1])/_l;
  d2f[1][1] = _h*(tmp2[0][2] + _stheta*tmp[2])/_l;
}

template <char p0, char p1, char p2>
void
F1_bv<p0, p1, p2, 2>::set_args_impl(double const u[3], double s_hat,
                                 double const s[3])
             
{
  _s_hat = s_hat;

  // _s[0] = s[0];
  // _s[1] = s[1];
  // _s[2] = s[2];
  _s0 = s[0];
  _ds[0] = s[1] - s[0];
  _ds[1] = s[2] - s[0];

  _u0 = u[0];
  _du[0] = u[1] - _u0;
  _du[1] = u[2] - _u0;
}

template <char p0, char p1, char p2>
void
F1_bv<p0, p1, p2, 2>::set_lambda_impl(double const lam[2])
{
  // double const lam0 = 1 - lam[0] - lam[1], lam1 = lam[0], lam2 = lam[1];

  // _stheta = (1 - _theta)*_s_hat + _theta*(lam0*_s[0] + lam1*_s[1] + lam2*_s[2]);
  _stheta = (1 - _theta)*_s_hat + _theta*(_s0 + _ds[0]*lam[0] + _ds[1]*lam[1]);

  _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];

  _p_lam[0] = component(p0, 0) + lam[0]*__dP(1, 0) + lam[1]*__dP(2, 0);
  _p_lam[1] = component(p0, 1) + lam[0]*__dP(1, 1) + lam[1]*__dP(2, 1);
  _p_lam[2] = component(p0, 2) + lam[0]*__dP(1, 2) + lam[1]*__dP(2, 2);

  _q = _p_lam[0]*_p_lam[0] + _p_lam[1]*_p_lam[1] + _p_lam[2]*_p_lam[2];

  _l = sqrt(_q);
}

#undef __dP
#undef __dP_dot_p_lam_n3

#endif // __COST_FUNCS_IMPL_HPP__
