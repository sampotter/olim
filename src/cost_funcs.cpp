#include "cost_funcs.hpp"

#include <cmath>

template <>
void
F0<3, 2>::eval_impl(double & f) const
{
  f = _u_lam + _sh*_l;
}

template <>
void
F0<3, 2>::grad_impl(double df[2]) const
{
  double const dPt_dot_p[2] = {
    _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2],
    _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2]
  };
  df[0] = _du[0] + _sh*dPt_dot_p[0]/_l;
  df[1] = _du[1] + _sh*dPt_dot_p[1]/_l;
}

template <>
void
F0<3, 2>::hess_impl(double d2f[3]) const
{
  double tmp[6]; // workspace
  tmp[0] = 1 - _p_lam[0]*_p_lam[0]/_q; // (0, 0)
  tmp[1] = -_p_lam[0]*_p_lam[1]/_q;    // (0, 1) & (1, 0)
  tmp[2] = -_p_lam[0]*_p_lam[2]/_q;    // (0, 2) & (2, 0)
  tmp[3] = 1 - _p_lam[1]*_p_lam[1]/_q; // (1, 1)
  tmp[4] = -_p_lam[1]*_p_lam[2]/_q;    // (1, 2) & (2, 1)
  tmp[5] = 1 - _p_lam[2]*_p_lam[2]/_q; // (2, 2)

  double tmp2[2][3];
  tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
  tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
  tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
  tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
  tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
  tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

  tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
  tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
  tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

  d2f[0] = _sh*tmp[0]/_l;
  d2f[1] = _sh*tmp[1]/_l;
  d2f[2] = _sh*tmp[2]/_l;
}

template <>
void
F0<3, 2>::set_lambda_impl(double const lam[2])
{
  _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];
  _p_lam[0] = _p0[0] + lam[0]*_dP[0][0] + lam[1]*_dP[1][0];
  _p_lam[1] = _p0[1] + lam[0]*_dP[0][1] + lam[1]*_dP[1][1];
  _p_lam[2] = _p0[2] + lam[0]*_dP[0][2] + lam[1]*_dP[1][2];
  _q = _p_lam[0]*_p_lam[0] + _p_lam[1]*_p_lam[1] + _p_lam[2]*_p_lam[2];
  _l = sqrt(_q);
}

template <>
void
F0<3, 2>::set_args_impl(double const u[3], double s_hat, double const s[3],
                        double const p[3][3])
{
  _sh = _h*((1 - _theta)*s_hat + _theta*(s[0] + s[1] + s[2])/3);
  _u0 = u[0];
  _du[0] = u[1] - _u0;
  _du[1] = u[2] - _u0;
  _p0[0] = p[0][0];
  _p0[1] = p[0][1];
  _p0[2] = p[0][2];
  _dP[0][0] = p[1][0] - p[0][0];
  _dP[0][1] = p[1][1] - p[0][1];
  _dP[0][2] = p[1][2] - p[0][2];
  _dP[1][0] = p[2][0] - p[0][0];
  _dP[1][1] = p[2][1] - p[0][1];
  _dP[1][2] = p[2][2] - p[0][2];
}

template <>
void
F1<3, 2>::eval_impl(double & f) const
{
  f = _u_lam + _h*_stheta*_l;
}

template <>
void
F1<3, 2>::grad_impl(double df[2]) const
{
  double const dPt_dot_p[2] = {
    _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2],
    _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2]
  };
  df[0] = _du[0] + _h*(_theta*_q*_ds[0] + _stheta*dPt_dot_p[0])/_l;
  df[1] = _du[1] + _h*(_theta*_q*_ds[1] + _stheta*dPt_dot_p[1])/_l;
}

template <>
void
F1<3, 2>::hess_impl(double d2f[3]) const
{
  double tmp[6]; // workspace
  tmp[0] = 1 - _p_lam[0]*_p_lam[0]/_q; // (0, 0)
  tmp[1] = -_p_lam[0]*_p_lam[1]/_q;    // (0, 1) & (1, 0)
  tmp[2] = -_p_lam[0]*_p_lam[2]/_q;    // (0, 2) & (2, 0)
  tmp[3] = 1 - _p_lam[1]*_p_lam[1]/_q; // (1, 1)
  tmp[4] = -_p_lam[1]*_p_lam[2]/_q;    // (1, 2) & (2, 1)
  tmp[5] = 1 - _p_lam[2]*_p_lam[2]/_q; // (2, 2)

  double tmp2[2][3];
  tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
  tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
  tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
  tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
  tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
  tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

  // compute dP'*(I - p*p'/q)*dP
  tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
  tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
  tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

  // compute dP'*p
  tmp[3] = _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2];
  tmp[4] = _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2];

  // compute theta*(dP'*p*ds' + ds*p'*dP)
  tmp2[0][0] = 2*_theta*tmp[3]*_ds[0];
  tmp2[0][1] = _theta*(tmp[3]*_ds[1] + tmp[4]*_ds[0]);
  tmp2[0][2] = 2*_theta*tmp[4]*_ds[1];

  // compute Hessian
  d2f[0] = _h*(tmp2[0][0] + _stheta*tmp[0])/_l;
  d2f[1] = _h*(tmp2[0][1] + _stheta*tmp[1])/_l;
  d2f[2] = _h*(tmp2[0][2] + _stheta*tmp[2])/_l;
}

template <>
void
F1<3, 2>::set_args_impl(double const u[3], double s_hat, double const s[3],
                        double const p[3][3])
{
  _s_hat = s_hat;
  _s[0] = s[0];
  _s[1] = s[1];
  _s[2] = s[2];
  _ds[0] = s[1] - s[0];
  _ds[1] = s[2] - s[0];
  _u0 = u[0];
  _du[0] = u[1] - _u0;
  _du[1] = u[2] - _u0;
  _p0[0] = p[0][0];
  _p0[1] = p[0][1];
  _p0[2] = p[0][2];
  _dP[0][0] = p[1][0] - p[0][0];
  _dP[0][1] = p[1][1] - p[0][1];
  _dP[0][2] = p[1][2] - p[0][2];
  _dP[1][0] = p[2][0] - p[0][0];
  _dP[1][1] = p[2][1] - p[0][1];
  _dP[1][2] = p[2][2] - p[0][2];
}

template <>
void
F1<3, 2>::set_lambda_impl(double const lam[2])
{
  double const lam0 = 1 - lam[0] - lam[1], lam1 = lam[0], lam2 = lam[1];
  _stheta = (1 - _theta)*_s_hat + _theta*(lam0*_s[0] + lam1*_s[1] + lam2*_s[2]);
  _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];
  _p_lam[0] = _p0[0] + lam[0]*_dP[0][0] + lam[1]*_dP[1][0];
  _p_lam[1] = _p0[1] + lam[0]*_dP[0][1] + lam[1]*_dP[1][1];
  _p_lam[2] = _p0[2] + lam[0]*_dP[0][2] + lam[1]*_dP[1][2];
  _q = _p_lam[0]*_p_lam[0] + _p_lam[1]*_p_lam[1] + _p_lam[2]*_p_lam[2];
  _l = sqrt(_q);
}

template <>
void
F0_fac<3, 2>::grad_impl(double df[2]) const
{
  double const dPt_dot_p[2] = {
    _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2],
    _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2]
  };
  double const dPt_dot_p_minus_p_fac[2] = {
    dPt_dot_p[0] - _dPt_dot_p_fac[0],
    dPt_dot_p[1] - _dPt_dot_p_fac[1]
  };
  df[0] = _dtau[0] + _sh*dPt_dot_p[0]/_l +
    _s_fac*_h*dPt_dot_p_minus_p_fac[0]/_l_fac_lam;
  df[1] = _dtau[1] + _sh*dPt_dot_p[1]/_l +
    _s_fac*_h*dPt_dot_p_minus_p_fac[1]/_l_fac_lam;
}

template <>
void
F0_fac<3, 2>::hess_impl(double d2f[3]) const
{
  F0<3, 2>::hess_impl(d2f);

  auto & z = _p_lam_minus_p_fac;
  auto & z_sq = _l_fac_lam_sq;

  double tmp[6]; // workspace
  tmp[0] = 1 - z[0]*z[0]/z_sq; // (0, 0)
  tmp[1] = -z[0]*z[1]/z_sq;    // (0, 1) & (1, 0)
  tmp[2] = -z[0]*z[2]/z_sq;    // (0, 2) & (2, 0)
  tmp[3] = 1 - z[1]*z[1]/z_sq; // (1, 1)
  tmp[4] = -z[1]*z[2]/z_sq;    // (1, 2) & (2, 1)
  tmp[5] = 1 - z[2]*z[2]/z_sq; // (2, 2)

  // TODO: the next two sections together are a good candidate for
  // being factored out as a function: they're just conjugating the
  // projector stored in tmp by dP (i.e. they compute dP'*P*dP in
  // place)

  double tmp2[2][3];
  tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
  tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
  tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
  tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
  tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
  tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

  tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
  tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
  tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

  d2f[0] += _s_fac*_h*tmp[0]/_l_fac_lam;
  d2f[1] += _s_fac*_h*tmp[1]/_l_fac_lam;
  d2f[2] += _s_fac*_h*tmp[2]/_l_fac_lam;
}

template <>
void
F0_fac<3, 2>::set_lambda_impl(double const lambda[2])
{
  F0<3, 2>::set_lambda_impl(lambda);

  _p_lam_minus_p_fac[0] = _p_fac[0] - _p_lam[0];
  _p_lam_minus_p_fac[1] = _p_fac[1] - _p_lam[1];
  _p_lam_minus_p_fac[2] = _p_fac[2] - _p_lam[2];

  auto const & z = _p_lam_minus_p_fac;

  _l_fac_lam_sq = z[0]*z[0] + z[1]*z[1] + z[2]*z[2];
  _l_fac_lam = std::sqrt(_l_fac_lam_sq);
}


template <>
void
F0_fac<3, 2>::set_args(double const u[3], double s_hat,
                       double const s[3], double const p[3][3],
                       double const p_fac[3], double s_fac)
{
  F0<3, 2>::set_args_impl(u, s_hat, s, p);

  memcpy(_p_fac, p_fac, 3*sizeof(double));

  _s_fac = s_fac;

  // These temporaries are for computing tau using T
  double l_fac[3], tmp[3];
  for (int i = 0; i < 3; ++i) {
    tmp[0] = p[i][0] - _p_fac[0];
    tmp[1] = p[i][1] - _p_fac[1];
    tmp[2] = p[i][2] - _p_fac[2];
    l_fac[i] = norm2<3>(tmp);
  }

  _tau0 = u[0] - _s_fac*_h*l_fac[0];
  _dtau[0] = u[1] - _s_fac*_h*l_fac[1];
  _dtau[1] = u[2] - _s_fac*_h*l_fac[2];

  _dPt_dot_p_fac[0] = _dP[0][0]*_p_fac[0] + _dP[0][1]*_p_fac[1] +
    _dP[0][2]*_p_fac[2];
  _dPt_dot_p_fac[1] = _dP[1][0]*_p_fac[0] + _dP[1][1]*_p_fac[1] +
    _dP[1][2]*_p_fac[2];
}

template <>
void
F1_fac<3, 2>::grad_impl(double df[2]) const
{
  double const dPt_dot_p[2] = {
    _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2],
    _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2]
  };
  double const dPt_dot_p_minus_p_fac[2] = {
    dPt_dot_p[0] - _dPt_dot_p_fac[0],
    dPt_dot_p[1] - _dPt_dot_p_fac[1]
  };
  df[0] = _dtau[0] + _h*(_theta*_q*_ds[0] + _stheta*dPt_dot_p[0])/_l +
    _s_fac*_h*dPt_dot_p_minus_p_fac[0]/_l_fac_lam;
  df[1] = _dtau[1] + _h*(_theta*_q*_ds[1] + _stheta*dPt_dot_p[1])/_l +
    _s_fac*_h*dPt_dot_p_minus_p_fac[1]/_l_fac_lam;
}

template <>
void
F1_fac<3, 2>::hess_impl(double d2f[3]) const
{
  F1<3, 2>::hess_impl(d2f);

  auto & z = _p_lam_minus_p_fac;
  auto & z_sq = _l_fac_lam_sq;

  double tmp[6]; // workspace
  tmp[0] = 1 - z[0]*z[0]/z_sq; // (0, 0)
  tmp[1] = -z[0]*z[1]/z_sq;    // (0, 1) & (1, 0)
  tmp[2] = -z[0]*z[2]/z_sq;    // (0, 2) & (2, 0)
  tmp[3] = 1 - z[1]*z[1]/z_sq; // (1, 1)
  tmp[4] = -z[1]*z[2]/z_sq;    // (1, 2) & (2, 1)
  tmp[5] = 1 - z[2]*z[2]/z_sq; // (2, 2)

  // TODO: the next two sections together are a good candidate for
  // being factored out as a function: they're just conjugating the
  // projector stored in tmp by dP (i.e. they compute dP'*P*dP in
  // place)

  double tmp2[2][3];
  tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
  tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
  tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
  tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
  tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
  tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

  tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
  tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
  tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

  d2f[0] += _s_fac*_h*tmp[0]/_l_fac_lam;
  d2f[1] += _s_fac*_h*tmp[1]/_l_fac_lam;
  d2f[2] += _s_fac*_h*tmp[2]/_l_fac_lam;
}

template <>
void
F1_fac<3, 2>::set_lambda_impl(double const lambda[2])
{
  F1<3, 2>::set_lambda_impl(lambda);

  _p_lam_minus_p_fac[0] = _p_fac[0] - _p_lam[0];
  _p_lam_minus_p_fac[1] = _p_fac[1] - _p_lam[1];
  _p_lam_minus_p_fac[2] = _p_fac[2] - _p_lam[2];

  auto const & z = _p_lam_minus_p_fac;

  _l_fac_lam_sq = z[0]*z[0] + z[1]*z[1] + z[2]*z[2];
  _l_fac_lam = std::sqrt(_l_fac_lam_sq);
}

template <>
void
F1_fac<3, 2>::set_args(double const u[3], double s_hat,
                       double const s[3], double const p[3][3],
                       double const p_fac[3], double s_fac)
{
  F1<3, 2>::set_args_impl(u, s_hat, s, p);

  memcpy(_p_fac, p_fac, 3*sizeof(double));

  _s_fac = s_fac;

  // These temporaries are for computing tau using T
  double l_fac[3], tmp[3];
  for (int i = 0; i < 3; ++i) {
    tmp[0] = p[i][0] - _p_fac[0];
    tmp[1] = p[i][1] - _p_fac[1];
    tmp[2] = p[i][2] - _p_fac[2];
    l_fac[i] = norm2<3>(tmp);
  }

  _tau0 = u[0] - _s_fac*_h*l_fac[0];
  _dtau[0] = u[1] - _s_fac*_h*l_fac[1];
  _dtau[1] = u[2] - _s_fac*_h*l_fac[2];

  _dPt_dot_p_fac[0] = _dP[0][0]*_p_fac[0] + _dP[0][1]*_p_fac[1] +
    _dP[0][2]*_p_fac[2];
  _dPt_dot_p_fac[1] = _dP[1][0]*_p_fac[0] + _dP[1][1]*_p_fac[1] +
    _dP[1][2]*_p_fac[2];
}
