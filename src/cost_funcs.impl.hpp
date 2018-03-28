#ifndef __COST_FUNCS_IMPL_HPP__
#define __COST_FUNCS_IMPL_HPP__

#include <cmath>

#include "update_rules.utils.hpp"

template <char pj, char p0>
constexpr char dp(char i) {
  return component(pj, i) - component(p0, i);
}

#define __dP(i, j) dp<p##i, p0>(j)

#define __dPt_dP(i, j) (                        \
    dp<p##i, p0>(0)*dp<p##j, p0>(0) +           \
    dp<p##i, p0>(1)*dp<p##j, p0>(1) +           \
    dp<p##i, p0>(2)*dp<p##j, p0>(2))

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
  df[0] = _du[0] + _sh*_y[0];
  df[1] = _du[1] + _sh*_y[1];
}

template <char p0, char p1, char p2>
void
F0_bv<p0, p1, p2, 2>::hess_impl(double d2f[3]) const
{
  d2f[0] = _sh*(__dPt_dP(1, 1) - _y[0]*_y[0])/_l;
  d2f[1] = _sh*(__dPt_dP(1, 2) - _y[0]*_y[1])/_l;
  d2f[2] = _sh*(__dPt_dP(2, 2) - _y[1]*_y[1])/_l;
}

template <char p0, char p1, char p2>
void
F0_bv<p0, p1, p2, 2>::set_lambda_impl(double const lam[2])
{
  _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];

  double const p[3] = {
    component(p0, 0) + lam[0]*__dP(1, 0) + lam[1]*__dP(2, 0),
    component(p0, 1) + lam[0]*__dP(1, 1) + lam[1]*__dP(2, 1),
    component(p0, 2) + lam[0]*__dP(1, 2) + lam[1]*__dP(2, 2)
  };
  _l = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

  _y[0] = (__dP(1, 0)*p[0] + __dP(1, 1)*p[1] + __dP(1, 2)*p[2])/_l;
  _y[1] = (__dP(2, 0)*p[0] + __dP(2, 1)*p[1] + __dP(2, 2)*p[2])/_l;
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
  f = _u_lam + _sh*_l;
}

template <char p0, char p1, char p2>
void
F1_bv<p0, p1, p2, 2>::grad_impl(double df[2]) const
{
  df[0] = _du[0] + _theta*_h*_l*_ds[0] + _sh*_y[0];
  df[1] = _du[1] + _theta*_h*_l*_ds[1] + _sh*_y[1];
}

template <char p0, char p1, char p2>
void
F1_bv<p0, p1, p2, 2>::hess_impl(double d2f[3]) const
{
  d2f[0] = _sh*(__dPt_dP(1, 1) - _y[0]*_y[0])/_l + 2*_h*_theta*_y[0]*_ds[0];
  d2f[1] = _sh*(__dPt_dP(1, 2) - _y[0]*_y[1])/_l + _h*_theta*(_y[0]*_ds[1] + _y[1]*_ds[0]);
  d2f[2] = _sh*(__dPt_dP(2, 2) - _y[1]*_y[1])/_l + 2*_h*_theta*_y[1]*_ds[1];
}

template <char p0, char p1, char p2>
void
F1_bv<p0, p1, p2, 2>::set_args_impl(double const u[3], double s_hat,
                                    double const s[3])
{
  _s_hat = s_hat;

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
  _sh = ((1 - _theta)*_s_hat + _theta*(_s0 + _ds[0]*lam[0] + _ds[1]*lam[1]))*_h;

  _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];

  double const p[3] = {
    component(p0, 0) + lam[0]*__dP(1, 0) + lam[1]*__dP(2, 0),
    component(p0, 1) + lam[0]*__dP(1, 1) + lam[1]*__dP(2, 1),
    component(p0, 2) + lam[0]*__dP(1, 2) + lam[1]*__dP(2, 2)
  };
  _l = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);

  _y[0] = (__dP(1, 0)*p[0] + __dP(1, 1)*p[1] + __dP(1, 2)*p[2])/_l;
  _y[1] = (__dP(2, 0)*p[0] + __dP(2, 1)*p[1] + __dP(2, 2)*p[2])/_l;
}

#undef __dP

#endif // __COST_FUNCS_IMPL_HPP__
