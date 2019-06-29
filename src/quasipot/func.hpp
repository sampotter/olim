#pragma once

#include "../mat.hpp"

template <int n, int d>
struct cost_functor
{
  vec<double, d> dU;
  mat<double, n, d> dX, dB;
  vec<double, n> x0, b0, xl, bl;
  double U0, Ul, xl_norm, bl_norm;

  vec<double, n> xl_over_xl_norm, bl_over_bl_norm;
  mat<double, n, d> dX_over_xl_norm, dB_over_bl_norm;

  cost_functor(vec<double, n> const * x, double const * U,
               vec<double, n> const * b) {
    x0 = x[0];
    U0 = U[0];
    b0 = b[0];

    for (int i = 0; i < d; ++i) {
      dX[i] = x[i + 1] - x0;
      dU[i] = U[i + 1] - U0;
      dB[i] = b[i + 1] - b0;
    }
  }

  inline void set_lambda(vec<double, d> const & lam) {
    Ul = U0 + dU*lam;

    xl = x0 + dX*lam;
    xl_norm = xl.norm2();
    xl_over_xl_norm = xl/xl_norm;
    dX_over_xl_norm = dX/xl_norm;

    bl = b0 + dB*lam;
    bl_norm = bl.norm2();
    bl_over_bl_norm = bl/bl_norm;
    dB_over_bl_norm = dB/bl_norm;
  }

  inline void eval(double & f) const {
    f = Ul + bl_norm*xl_norm - bl*xl;
  }

  // TODO: below, we can *DEFINITELY* cache a few more things, but
  // let's just try to get this working, first

  inline void grad(vec<double, d> & df) const {
    df = dU + bl_norm*xl_norm*(dB_over_bl_norm - dX_over_xl_norm).t()*(
      bl_over_bl_norm - xl_over_xl_norm);
  }

  // TODO: want to add a "mat" class and in particular a "symmat"
  // class...
  //
  // - we will also want to make "sym_inner" and "sym_outer" functions which
  //   take matrices and return symmetric
  inline void hess(mat<double, d, d> & d2f) const {
    auto tmp1 = dB_over_bl_norm - dX_over_xl_norm;
    auto tmp2 = dB_over_bl_norm*bl_over_bl_norm - dX_over_xl_norm*xl_over_xl_norm;
    d2f = bl_norm*xl_norm*(tmp1.t()*tmp1 - tmp2*tmp2.t());
  }
};
