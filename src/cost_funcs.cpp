#include "cost_funcs.hpp"

#include <cassert>
#include <cmath>
#include <utility>

#define __invert2x2inplace(X) do {              \
    double const det = X[0]*X[2] - X[1]*X[1];   \
    std::swap(X[0], X[2]);                      \
    X[0] /= det;                                \
    X[1] /= -det;                               \
    X[2] /= det;                                \
  } while (0)

template <>
void lagmults<2>(double const * lam, double const * df, double const * d2f,
                 double * mu, int * k)
{
  // First, find the active set for lam

  int active_set[2] = {-1, -1};
  int num_active = 0;

  if (fabs(lam[1]) < EPS(double)) active_set[num_active++] = 1;
  if (fabs(lam[0]) < EPS(double)) active_set[num_active++] = 0;
  if (fabs(1 - lam[0] - lam[1]) < EPS(double)) active_set[num_active++] = 2;

  *k = num_active;

  // If the active_set is {0, 1}, then mu = -df, which means that we
  // can return early.

  if (num_active == 2 &&
      ((active_set[0] == 0 && active_set[1] == 1) ||
       (active_set[0] == 1 && active_set[1] == 0))) {
    mu[0] = -df[0];
    mu[1] = -df[1];
  } else {
    double inv_d2f[3] = {d2f[0], d2f[1], d2f[2]};
    __invert2x2inplace(inv_d2f);

    double p[2];
    p[0] = inv_d2f[0]*df[0] + inv_d2f[1]*df[1];
    p[1] = inv_d2f[1]*df[0] + inv_d2f[2]*df[1];

    if (num_active == 1) {
      if (active_set[0] == 0) {
        mu[0] = -p[0]/d2f[0];
      } else if (active_set[0] == 1) {
        mu[0] = -p[1]/d2f[2];
      } else if (active_set[0] == 2) {
        mu[0] = (p[0] + p[1])/(d2f[0] + 2*d2f[1] + d2f[2]);
      } else {
        assert(false);
        mu[0] = mu[1] = INF(double);
      }
    } else if (num_active == 2) {
      // TODO: we can do this in place without defining A and b
      double A[3], b[2];
      if ((active_set[0] == 0 && active_set[1] == 2) ||
          (active_set[0] == 2 && active_set[1] == 0)) {
        A[0] = d2f[0];
        A[1] = -(d2f[0] + d2f[1]);
        A[2] = d2f[0] + 2*d2f[1] + d2f[2];
        __invert2x2inplace(A);
        b[0] = -p[0];
        b[1] = p[0] + p[1];
      } else if ((active_set[0] == 1 && active_set[1] == 2) ||
                 (active_set[0] == 2 && active_set[1] == 1)) {
        A[0] = d2f[2];
        A[1] = -(d2f[1] + d2f[2]);
        A[2] = d2f[0] + 2*d2f[1] + d2f[2];
        __invert2x2inplace(A);
        b[0] = -p[1];
        b[1] = p[0] + p[1];
      } else {
        assert(false);
        A[0] = A[1] = A[2] = b[0] = b[1] = INF(double);
      }
      mu[0] = A[0]*b[0] + A[1]*b[1];
      mu[1] = A[1]*b[0] + A[2]*b[1];
    } else {
      assert(false);
      mu[0] = mu[1] = INF(double);
    }
  }

  __check(mu[0]);
  __check(mu[1]);
}

#undef __invert2x2inplace


void eval(eval_wkspc const & w, double & f)
{
  f = w.u_lam + w.sh*w.l_lam;
  __check(f);
}

void eval(eval_wkspc const & w, fac_eval_wkspc const & fw, double & f)
{
  f = fw.tau_lam + fw.sh_fac*fw.l_fac_lam + w.sh*w.l_lam;
  __check(f);
}

// template <>
// void grad<F0_wkspc<3, 2>>(F0_wkspc<3, 2> const & w, double * df)
// {
//   df[0] = w.du[0] + w.sh*w.dPt_nu_lam[0];
//   df[1] = w.du[1] + w.sh*w.dPt_nu_lam[1];
//   __check(df[0]);
//   __check(df[1]);
// }

// template <>
// void grad<F0_wkspc<3, 2>>(
//   F0_wkspc<3, 2> const & w, fac_wkspc const & fw, double * df)
// {
//   df[0] = fw.dtau[0] + w.sh*w.dPt_nu_lam[0];
//   df[1] = fw.dtau[0] + w.sh*w.dPt_nu_lam[1];

//   if (fw.l_fac_lam > 1e1*EPS(double)) {
//     df[0] += fw.sh_fac*fw.dPt_nu_fac_lam[0];
//     df[1] += fw.sh_fac*fw.dPt_nu_fac_lam[1];
//   }

//   __check(df[0]);
//   __check(df[1]);
// }

// template <>
// void grad<F1_wkspc<3, 2>>(F1_wkspc<3, 2> const & w, double * df)
// {
//   df[0] = w.du[0] + w.theta_h_ds[0]*w.l_lam + w.sh*w.dPt_nu_lam[0];
//   df[1] = w.du[1] + w.theta_h_ds[1]*w.l_lam + w.sh*w.dPt_nu_lam[1];
//   __check(df[0]);
//   __check(df[1]);
// }

// template <>
// void grad<F1_wkspc<3, 2>>(F1_wkspc<3, 2> const & w, double * df)
// {
//   df[0] = fw.dtau[0] + w.theta_h_ds[0]*w.l_lam + w.sh*w.dPt_nu_lam[0];
//   df[1] = fw.dtau[1] + w.theta_h_ds[1]*w.l_lam + w.sh*w.dPt_nu_lam[1];

//   if (fw.l_fac_lam > 1e1*EPS(double)) {
//     df[0] += fw.sh_fac*fw.dPt_nu_fac_lam[0];
//     df[1] += fw.sh_fac*fw.dPt_nu_fac_lam[1];
//   }

//   __check(df[0]);
//   __check(df[1]);
// }

// template <>
// void hess<F0_wkspc<3, 2>>(F0_wkspc<3, 2> const & w, double * d2f)
// {
//   double tmp = w.sh/w.l_lam;

//   d2f[0] = tmp*(w.dPt_dP[0] - w.dPt_nu_lam[0]*w.dPt_nu_lam[0]);
//   d2f[1] = tmp*(w.dPt_dP[1] - w.dPt_nu_lam[0]*w.dPt_nu_lam[1]);
//   d2f[2] = tmp*(w.dPt_dP[2] - w.dPt_nu_lam[1]*w.dPt_nu_lam[1]);

//   __check(d2f[0]);
//   __check(d2f[1]);
//   __check(d2f[2]);
// }

// template <>
// void hess<F0_wkspc<3, 2>>(
//   F0_wkspc<3, 2> const & w, fac_wkspc const & fw, double * d2f)
// {
//   hess(w, d2f);

//   if (_l_fac_lam > EPS(double)) {
//     double tmp = fw.sh_fac/fw.l_fac_lam;
//     d2f[0] += tmp*(w.dPt_dP[0] - w.dPt_nu_fac_lam[0]*w.dPt_nu_fac_lam[0]);
//     d2f[1] += tmp*(w.dPt_dP[1] - w.dPt_nu_fac_lam[0]*w.dPt_nu_fac_lam[1]);
//     d2f[2] += tmp*(w.dPt_dP[2] - w.dPt_nu_fac_lam[1]*w.dPt_nu_fac_lam[1]);
//   }

//   __check(d2f[0]);
//   __check(d2f[1]);
//   __check(d2f[2]);
// }

// template <>
// void hess<F1_wkspc<3, 2>>(F1_wkspc<3, 2> const & w, double * d2f)
// {
//   double tmp = w.sh/w.l_lam;

//   d2f[0] = 2*w.theta_h_ds[0]*w.dPt_nu_lam[0] +
//     tmp*(w.dPt_dP[0] - w.dPt_nu_lam[0]*w.dPt_nu_lam[0]);
//   d2f[1] = w.theta_h_ds[0]*w.dPt_nu_lam[1] + w.theta_h_ds[1]*w.dPt_nu_lam[0]
//     tmp*(w.dPt_dP[1] - w.dPt_nu_lam[0]*w.dPt_nu_lam[1]);
//   d2f[2] = 2*w.theta_h_ds[1]*w.dPt_nu_lam[1] +
//     tmp*(w.dPt_dP[2] - w.dPt_nu_lam[1]*w.dPt_nu_lam[1]);

//   __check(d2f[0]);
//   __check(d2f[1]);
//   __check(d2f[2]);
// }

// template <>
// void hess<F1_wkspc<3, 2>>(
//   F1_wkspc<3, 2> const & w, fac_wkspc const & fw, double * d2f)
// {
//   hess(w, d2f);

//   if (_l_fac_lam > EPS(double)) {
//     double tmp = fw.sh_fac/fw.l_fac_lam;
//     d2f[0] += tmp*(w.dPt_dP[0] - w.dPt_nu_fac_lam[0]*w.dPt_nu_fac_lam[0]);
//     d2f[1] += tmp*(w.dPt_dP[1] - w.dPt_nu_fac_lam[0]*w.dPt_nu_fac_lam[1]);
//     d2f[2] += tmp*(w.dPt_dP[2] - w.dPt_nu_fac_lam[1]*w.dPt_nu_fac_lam[1]);
//   }

//   __check(d2f[0]);
//   __check(d2f[1]);
//   __check(d2f[2]);
// }

////////////////////////////////////////////////////////////////////////////////
// OLD /////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// template <>
// void
// F0<3, 2>::eval_impl(double & f) const
// {
//   f = _u_lam + _sh*_l;
// }

// template <>
// void
// F0<3, 2>::grad_impl(double df[2]) const
// {
//   double const dPt_dot_p[2] = {
//     _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2],
//     _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2]
//   };
//   df[0] = _du[0] + _sh*dPt_dot_p[0]/_l;
//   df[1] = _du[1] + _sh*dPt_dot_p[1]/_l;
// }

// template <>
// void
// F0<3, 2>::hess_impl(double d2f[3]) const
// {
//   double tmp[6]; // workspace
//   tmp[0] = 1 - _p_lam[0]*_p_lam[0]/_q; // (0, 0)
//   tmp[1] = -_p_lam[0]*_p_lam[1]/_q;    // (0, 1) & (1, 0)
//   tmp[2] = -_p_lam[0]*_p_lam[2]/_q;    // (0, 2) & (2, 0)
//   tmp[3] = 1 - _p_lam[1]*_p_lam[1]/_q; // (1, 1)
//   tmp[4] = -_p_lam[1]*_p_lam[2]/_q;    // (1, 2) & (2, 1)
//   tmp[5] = 1 - _p_lam[2]*_p_lam[2]/_q; // (2, 2)

//   double tmp2[2][3];
//   tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
//   tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
//   tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
//   tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
//   tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
//   tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

//   tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
//   tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
//   tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

//   d2f[0] = _sh*tmp[0]/_l;
//   d2f[1] = _sh*tmp[1]/_l;
//   d2f[2] = _sh*tmp[2]/_l;
// }

// template <>
// void
// F0<3, 2>::set_lambda_impl(double const lam[2])
// {
//   _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];
//   _p_lam[0] = _p0[0] + lam[0]*_dP[0][0] + lam[1]*_dP[1][0];
//   _p_lam[1] = _p0[1] + lam[0]*_dP[0][1] + lam[1]*_dP[1][1];
//   _p_lam[2] = _p0[2] + lam[0]*_dP[0][2] + lam[1]*_dP[1][2];
//   _q = _p_lam[0]*_p_lam[0] + _p_lam[1]*_p_lam[1] + _p_lam[2]*_p_lam[2];
//   _l = sqrt(_q);
// }

// template <>
// void
// F0<3, 2>::set_args_impl(double const u[3], double s_hat, double const s[3],
//                         double const p[3][3])
// {
//   _sh = _h*((1 - _theta)*s_hat + _theta*(s[0] + s[1] + s[2])/3);
//   _u0 = u[0];
//   _du[0] = u[1] - _u0;
//   _du[1] = u[2] - _u0;
//   _p0[0] = p[0][0];
//   _p0[1] = p[0][1];
//   _p0[2] = p[0][2];
//   _dP[0][0] = p[1][0] - p[0][0];
//   _dP[0][1] = p[1][1] - p[0][1];
//   _dP[0][2] = p[1][2] - p[0][2];
//   _dP[1][0] = p[2][0] - p[0][0];
//   _dP[1][1] = p[2][1] - p[0][1];
//   _dP[1][2] = p[2][2] - p[0][2];
// }

// template <>
// void
// F1<3, 2>::eval_impl(double & f) const
// {
//   f = _u_lam + _h*_stheta*_l;
// }

// template <>
// void
// F1<3, 2>::grad_impl(double df[2]) const
// {
//   double const dPt_dot_p[2] = {
//     _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2],
//     _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2]
//   };
//   df[0] = _du[0] + _h*(_theta*_q*_ds[0] + _stheta*dPt_dot_p[0])/_l;
//   df[1] = _du[1] + _h*(_theta*_q*_ds[1] + _stheta*dPt_dot_p[1])/_l;
// }

// template <>
// void
// F1<3, 2>::hess_impl(double d2f[3]) const
// {
//   double tmp[6]; // workspace
//   tmp[0] = 1 - _p_lam[0]*_p_lam[0]/_q; // (0, 0)
//   tmp[1] = -_p_lam[0]*_p_lam[1]/_q;    // (0, 1) & (1, 0)
//   tmp[2] = -_p_lam[0]*_p_lam[2]/_q;    // (0, 2) & (2, 0)
//   tmp[3] = 1 - _p_lam[1]*_p_lam[1]/_q; // (1, 1)
//   tmp[4] = -_p_lam[1]*_p_lam[2]/_q;    // (1, 2) & (2, 1)
//   tmp[5] = 1 - _p_lam[2]*_p_lam[2]/_q; // (2, 2)

//   double tmp2[2][3];
//   tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
//   tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
//   tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
//   tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
//   tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
//   tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

//   // compute dP'*(I - p*p'/q)*dP
//   tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
//   tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
//   tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

//   // compute dP'*p
//   tmp[3] = _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2];
//   tmp[4] = _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2];

//   // compute theta*(dP'*p*ds' + ds*p'*dP)
//   tmp2[0][0] = 2*_theta*tmp[3]*_ds[0];
//   tmp2[0][1] = _theta*(tmp[3]*_ds[1] + tmp[4]*_ds[0]);
//   tmp2[0][2] = 2*_theta*tmp[4]*_ds[1];

//   // compute Hessian
//   d2f[0] = _h*(tmp2[0][0] + _stheta*tmp[0])/_l;
//   d2f[1] = _h*(tmp2[0][1] + _stheta*tmp[1])/_l;
//   d2f[2] = _h*(tmp2[0][2] + _stheta*tmp[2])/_l;
// }

// template <>
// void
// F1<3, 2>::set_args_impl(double const u[3], double s_hat, double const s[3],
//                         double const p[3][3])
// {
//   _s_hat = s_hat;
//   _s[0] = s[0];
//   _s[1] = s[1];
//   _s[2] = s[2];
//   _ds[0] = s[1] - s[0];
//   _ds[1] = s[2] - s[0];
//   _u0 = u[0];
//   _du[0] = u[1] - _u0;
//   _du[1] = u[2] - _u0;
//   _p0[0] = p[0][0];
//   _p0[1] = p[0][1];
//   _p0[2] = p[0][2];
//   _dP[0][0] = p[1][0] - p[0][0];
//   _dP[0][1] = p[1][1] - p[0][1];
//   _dP[0][2] = p[1][2] - p[0][2];
//   _dP[1][0] = p[2][0] - p[0][0];
//   _dP[1][1] = p[2][1] - p[0][1];
//   _dP[1][2] = p[2][2] - p[0][2];
// }

// template <>
// void
// F1<3, 2>::set_lambda_impl(double const lam[2])
// {
//   double const lam0 = 1 - lam[0] - lam[1], lam1 = lam[0], lam2 = lam[1];
//   _stheta = (1 - _theta)*_s_hat + _theta*(lam0*_s[0] + lam1*_s[1] + lam2*_s[2]);
//   _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];
//   _p_lam[0] = _p0[0] + lam[0]*_dP[0][0] + lam[1]*_dP[1][0];
//   _p_lam[1] = _p0[1] + lam[0]*_dP[0][1] + lam[1]*_dP[1][1];
//   _p_lam[2] = _p0[2] + lam[0]*_dP[0][2] + lam[1]*_dP[1][2];
//   _q = _p_lam[0]*_p_lam[0] + _p_lam[1]*_p_lam[1] + _p_lam[2]*_p_lam[2];
//   _l = sqrt(_q);
// }

// template <>
// void
// F0_fac<3, 2>::eval_impl(double & f) const
// {
//   f = _tau_lam + _s_fac*_h*_l_fac_lam + _sh*_l;
// }

// template <>
// void
// F0_fac<3, 2>::grad_impl(double df[2]) const
// {
//   double const dPt_dot_p[2] = {
//     _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2],
//     _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2]
//   };

//   double const dPt_dot_p_minus_p_fac[2] = {
//     dPt_dot_p[0] - _dPt_dot_p_fac[0],
//     dPt_dot_p[1] - _dPt_dot_p_fac[1]
//   };

//   df[0] = _dtau[0] + _sh*dPt_dot_p[0]/_l;
//   df[1] = _dtau[1] + _sh*dPt_dot_p[1]/_l;

//   if (_l_fac_lam > 1e1*EPS(double)) {
//     df[0] += _s_fac*_h*dPt_dot_p_minus_p_fac[0]/_l_fac_lam;
//     df[1] += _s_fac*_h*dPt_dot_p_minus_p_fac[1]/_l_fac_lam;
//   }

//   __check(df[0]);
//   __check(df[1]);
// }

// template <>
// void
// F0_fac<3, 2>::hess_impl(double d2f[3]) const
// {
//   {
//     double tmp[6]; // workspace
//     tmp[0] = 1 - _p_lam[0]*_p_lam[0]/_q; // (0, 0)
//     tmp[1] = -_p_lam[0]*_p_lam[1]/_q;    // (0, 1) & (1, 0)
//     tmp[2] = -_p_lam[0]*_p_lam[2]/_q;    // (0, 2) & (2, 0)
//     tmp[3] = 1 - _p_lam[1]*_p_lam[1]/_q; // (1, 1)
//     tmp[4] = -_p_lam[1]*_p_lam[2]/_q;    // (1, 2) & (2, 1)
//     tmp[5] = 1 - _p_lam[2]*_p_lam[2]/_q; // (2, 2)

//     double tmp2[2][3];
//     tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
//     tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
//     tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
//     tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
//     tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
//     tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

//     tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
//     tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
//     tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

//     d2f[0] = _sh*tmp[0]/_l;
//     d2f[1] = _sh*tmp[1]/_l;
//     d2f[2] = _sh*tmp[2]/_l;
//   }

//   if (_l_fac_lam > EPS(double)) {
//     auto & z = _p_lam_minus_p_fac;
//     auto & z_sq = _l_fac_lam_sq;

//     double tmp[6]; // workspace
//     tmp[0] = 1 - z[0]*z[0]/z_sq; // (0, 0)
//     tmp[1] = -z[0]*z[1]/z_sq;    // (0, 1) & (1, 0)
//     tmp[2] = -z[0]*z[2]/z_sq;    // (0, 2) & (2, 0)
//     tmp[3] = 1 - z[1]*z[1]/z_sq; // (1, 1)
//     tmp[4] = -z[1]*z[2]/z_sq;    // (1, 2) & (2, 1)
//     tmp[5] = 1 - z[2]*z[2]/z_sq; // (2, 2)

//     // TODO: the next two sections together are a good candidate for
//     // being factored out as a function: they're just conjugating the
//     // projector stored in tmp by dP (i.e. they compute dP'*P*dP in
//     // place)

//     double tmp2[2][3];
//     tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
//     tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
//     tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
//     tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
//     tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
//     tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

//     tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
//     tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
//     tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

//     d2f[0] += _s_fac*_h*tmp[0]/_l_fac_lam;
//     d2f[1] += _s_fac*_h*tmp[1]/_l_fac_lam;
//     d2f[2] += _s_fac*_h*tmp[2]/_l_fac_lam;
//   }
// }

// template <>
// void
// F0_fac<3, 2>::set_lambda_impl(double const lam[2])
// {
//   _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];
//   _p_lam[0] = _p0[0] + lam[0]*_dP[0][0] + lam[1]*_dP[1][0];
//   _p_lam[1] = _p0[1] + lam[0]*_dP[0][1] + lam[1]*_dP[1][1];
//   _p_lam[2] = _p0[2] + lam[0]*_dP[0][2] + lam[1]*_dP[1][2];
//   _q = _p_lam[0]*_p_lam[0] + _p_lam[1]*_p_lam[1] + _p_lam[2]*_p_lam[2];
//   _l = sqrt(_q);

//   _p_lam_minus_p_fac[0] = _p_fac[0] - _p_lam[0];
//   _p_lam_minus_p_fac[1] = _p_fac[1] - _p_lam[1];
//   _p_lam_minus_p_fac[2] = _p_fac[2] - _p_lam[2];

//   auto const & z = _p_lam_minus_p_fac;

//   _l_fac_lam_sq = z[0]*z[0] + z[1]*z[1] + z[2]*z[2];
//   _l_fac_lam = std::sqrt(_l_fac_lam_sq);

//   _tau_lam = _tau0 + _dtau[0]*lam[0] + _dtau[1]*lam[1];
// }


// template <>
// void
// F0_fac<3, 2>::set_args(double const u[3], double s_hat,
//                        double const s[3], double const p[3][3],
//                        double const p_fac[3], double s_fac)
// {
//   _sh = _h*((1 - _theta)*s_hat + _theta*(s[0] + s[1] + s[2])/3);
//   _u0 = u[0];
//   _du[0] = u[1] - _u0;
//   _du[1] = u[2] - _u0;
//   _p0[0] = p[0][0];
//   _p0[1] = p[0][1];
//   _p0[2] = p[0][2];
//   _dP[0][0] = p[1][0] - p[0][0];
//   _dP[0][1] = p[1][1] - p[0][1];
//   _dP[0][2] = p[1][2] - p[0][2];
//   _dP[1][0] = p[2][0] - p[0][0];
//   _dP[1][1] = p[2][1] - p[0][1];
//   _dP[1][2] = p[2][2] - p[0][2];

//   memcpy(_p_fac, p_fac, 3*sizeof(double));

//   _s_fac = s_fac;

//   // These temporaries are for computing tau using T
//   double l_fac[3], tmp[3];
//   for (int i = 0; i < 3; ++i) {
//     tmp[0] = p[i][0] - _p_fac[0];
//     tmp[1] = p[i][1] - _p_fac[1];
//     tmp[2] = p[i][2] - _p_fac[2];
//     l_fac[i] = norm2<3>(tmp);
//   }

//   _tau0 = u[0] - _s_fac*_h*l_fac[0];
//   _dtau[0] = u[1] - _s_fac*_h*l_fac[1] - _tau0;
//   _dtau[1] = u[2] - _s_fac*_h*l_fac[2] - _tau0;

//   _dPt_dot_p_fac[0] = _dP[0][0]*_p_fac[0] + _dP[0][1]*_p_fac[1] +
//     _dP[0][2]*_p_fac[2];
//   _dPt_dot_p_fac[1] = _dP[1][0]*_p_fac[0] + _dP[1][1]*_p_fac[1] +
//     _dP[1][2]*_p_fac[2];
// }

// template <>
// void
// F1_fac<3, 2>::eval_impl(double & f) const
// {
//   f = _tau_lam + _s_fac*_h*_l_fac_lam + _stheta*_h*_l;
// }

// template <>
// void
// F1_fac<3, 2>::grad_impl(double df[2]) const
// {
//   double const dPt_dot_p[2] = {
//     _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2],
//     _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2]
//   };

//   double const dPt_dot_p_minus_p_fac[2] = {
//     dPt_dot_p[0] - _dPt_dot_p_fac[0],
//     dPt_dot_p[1] - _dPt_dot_p_fac[1]
//   };

//   df[0] = _dtau[0] + _h*(_theta*_q*_ds[0] + _stheta*dPt_dot_p[0])/_l;
//   df[1] = _dtau[1] + _h*(_theta*_q*_ds[1] + _stheta*dPt_dot_p[1])/_l;

//   if (_l_fac_lam > EPS(double)) {
//     df[0] += _s_fac*_h*dPt_dot_p_minus_p_fac[0]/_l_fac_lam;
//     df[1] += _s_fac*_h*dPt_dot_p_minus_p_fac[1]/_l_fac_lam;
//   }
// }

// template <>
// void
// F1_fac<3, 2>::hess_impl(double d2f[3]) const
// {
//   {
//     double tmp[6]; // workspace
//     tmp[0] = 1 - _p_lam[0]*_p_lam[0]/_q; // (0, 0)
//     tmp[1] = -_p_lam[0]*_p_lam[1]/_q;    // (0, 1) & (1, 0)
//     tmp[2] = -_p_lam[0]*_p_lam[2]/_q;    // (0, 2) & (2, 0)
//     tmp[3] = 1 - _p_lam[1]*_p_lam[1]/_q; // (1, 1)
//     tmp[4] = -_p_lam[1]*_p_lam[2]/_q;    // (1, 2) & (2, 1)
//     tmp[5] = 1 - _p_lam[2]*_p_lam[2]/_q; // (2, 2)

//     double tmp2[2][3];
//     tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
//     tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
//     tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
//     tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
//     tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
//     tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

//     // compute dP'*(I - p*p'/q)*dP
//     tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
//     tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
//     tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

//     // compute dP'*p
//     tmp[3] = _dP[0][0]*_p_lam[0] + _dP[0][1]*_p_lam[1] + _dP[0][2]*_p_lam[2];
//     tmp[4] = _dP[1][0]*_p_lam[0] + _dP[1][1]*_p_lam[1] + _dP[1][2]*_p_lam[2];

//     // compute theta*(dP'*p*ds' + ds*p'*dP)
//     tmp2[0][0] = 2*_theta*tmp[3]*_ds[0];
//     tmp2[0][1] = _theta*(tmp[3]*_ds[1] + tmp[4]*_ds[0]);
//     tmp2[0][2] = 2*_theta*tmp[4]*_ds[1];

//     // compute Hessian
//     d2f[0] = _h*(tmp2[0][0] + _stheta*tmp[0])/_l;
//     d2f[1] = _h*(tmp2[0][1] + _stheta*tmp[1])/_l;
//     d2f[2] = _h*(tmp2[0][2] + _stheta*tmp[2])/_l;
//   }

//   if (_l_fac_lam > EPS(double)) {
//     auto & z = _p_lam_minus_p_fac;
//     auto & z_sq = _l_fac_lam_sq;

//     double tmp[6]; // workspace
//     tmp[0] = 1 - z[0]*z[0]/z_sq; // (0, 0)
//     tmp[1] = -z[0]*z[1]/z_sq;    // (0, 1) & (1, 0)
//     tmp[2] = -z[0]*z[2]/z_sq;    // (0, 2) & (2, 0)
//     tmp[3] = 1 - z[1]*z[1]/z_sq; // (1, 1)
//     tmp[4] = -z[1]*z[2]/z_sq;    // (1, 2) & (2, 1)
//     tmp[5] = 1 - z[2]*z[2]/z_sq; // (2, 2)

//     // TODO: the next two sections together are a good candidate for
//     // being factored out as a function: they're just conjugating the
//     // projector stored in tmp by dP (i.e. they compute dP'*P*dP in
//     // place)

//     double tmp2[2][3];
//     tmp2[0][0] = tmp[0]*_dP[0][0] + tmp[1]*_dP[0][1] + tmp[2]*_dP[0][2];
//     tmp2[0][1] = tmp[1]*_dP[0][0] + tmp[3]*_dP[0][1] + tmp[4]*_dP[0][2];
//     tmp2[0][2] = tmp[2]*_dP[0][0] + tmp[4]*_dP[0][1] + tmp[5]*_dP[0][2];
//     tmp2[1][0] = tmp[0]*_dP[1][0] + tmp[1]*_dP[1][1] + tmp[2]*_dP[1][2];
//     tmp2[1][1] = tmp[1]*_dP[1][0] + tmp[3]*_dP[1][1] + tmp[4]*_dP[1][2];
//     tmp2[1][2] = tmp[2]*_dP[1][0] + tmp[4]*_dP[1][1] + tmp[5]*_dP[1][2];

//     tmp[0] = _dP[0][0]*tmp2[0][0] + _dP[0][1]*tmp2[0][1] + _dP[0][2]*tmp2[0][2];
//     tmp[1] = _dP[0][0]*tmp2[1][0] + _dP[0][1]*tmp2[1][1] + _dP[0][2]*tmp2[1][2];
//     tmp[2] = _dP[1][0]*tmp2[1][0] + _dP[1][1]*tmp2[1][1] + _dP[1][2]*tmp2[1][2];

//     d2f[0] += _s_fac*_h*tmp[0]/_l_fac_lam;
//     d2f[1] += _s_fac*_h*tmp[1]/_l_fac_lam;
//     d2f[2] += _s_fac*_h*tmp[2]/_l_fac_lam;
//   }
// }

// template <>
// void
// F1_fac<3, 2>::set_lambda_impl(double const lam[2])
// {
//   double const lam0 = 1 - lam[0] - lam[1], lam1 = lam[0], lam2 = lam[1];
//   _stheta = (1 - _theta)*_s_hat + _theta*(lam0*_s[0] + lam1*_s[1] + lam2*_s[2]);
//   _u_lam = _u0 + lam[0]*_du[0] + lam[1]*_du[1];
//   _p_lam[0] = _p0[0] + lam[0]*_dP[0][0] + lam[1]*_dP[1][0];
//   _p_lam[1] = _p0[1] + lam[0]*_dP[0][1] + lam[1]*_dP[1][1];
//   _p_lam[2] = _p0[2] + lam[0]*_dP[0][2] + lam[1]*_dP[1][2];
//   _q = _p_lam[0]*_p_lam[0] + _p_lam[1]*_p_lam[1] + _p_lam[2]*_p_lam[2];
//   _l = sqrt(_q);

//   _p_lam_minus_p_fac[0] = _p_fac[0] - _p_lam[0];
//   _p_lam_minus_p_fac[1] = _p_fac[1] - _p_lam[1];
//   _p_lam_minus_p_fac[2] = _p_fac[2] - _p_lam[2];

//   auto const & z = _p_lam_minus_p_fac;

//   _l_fac_lam_sq = z[0]*z[0] + z[1]*z[1] + z[2]*z[2];
//   _l_fac_lam = std::sqrt(_l_fac_lam_sq);

//   _tau_lam = _tau0 + _dtau[0]*lam[0] + _dtau[1]*lam[1];
// }

// template <>
// void
// F1_fac<3, 2>::set_args(double const u[3], double s_hat,
//                        double const s[3], double const p[3][3],
//                        double const p_fac[3], double s_fac)
// {
//   _s_hat = s_hat;
//   _s[0] = s[0];
//   _s[1] = s[1];
//   _s[2] = s[2];
//   _ds[0] = s[1] - s[0];
//   _ds[1] = s[2] - s[0];
//   _u0 = u[0];
//   _du[0] = u[1] - _u0;
//   _du[1] = u[2] - _u0;
//   _p0[0] = p[0][0];
//   _p0[1] = p[0][1];
//   _p0[2] = p[0][2];
//   _dP[0][0] = p[1][0] - p[0][0];
//   _dP[0][1] = p[1][1] - p[0][1];
//   _dP[0][2] = p[1][2] - p[0][2];
//   _dP[1][0] = p[2][0] - p[0][0];
//   _dP[1][1] = p[2][1] - p[0][1];
//   _dP[1][2] = p[2][2] - p[0][2];

//   memcpy(_p_fac, p_fac, 3*sizeof(double));

//   _s_fac = s_fac;

//   // These temporaries are for computing tau using T
//   double l_fac[3], tmp[3];
//   for (int i = 0; i < 3; ++i) {
//     tmp[0] = p[i][0] - _p_fac[0];
//     tmp[1] = p[i][1] - _p_fac[1];
//     tmp[2] = p[i][2] - _p_fac[2];
//     l_fac[i] = norm2<3>(tmp);
//   }

//   _tau0 = u[0] - _s_fac*_h*l_fac[0];
//   _dtau[0] = u[1] - _s_fac*_h*l_fac[1] - _tau0;
//   _dtau[1] = u[2] - _s_fac*_h*l_fac[2] - _tau0;

//   _dPt_dot_p_fac[0] = _dP[0][0]*_p_fac[0] + _dP[0][1]*_p_fac[1] +
//     _dP[0][2]*_p_fac[2];
//   _dPt_dot_p_fac[1] = _dP[1][0]*_p_fac[0] + _dP[1][1]*_p_fac[1] +
//     _dP[1][2]*_p_fac[2];
// }
