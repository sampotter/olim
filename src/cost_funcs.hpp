#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

#include <src/config.hpp>

#include "common.macros.hpp"

#include <cassert>
#include <cmath>

#define __check(x) do {                         \
    assert(!std::isinf(x));                     \
    assert(!std::isnan(x));                     \
  } while (0)

#define __sym_mat_size(d) ((d*(d + 1))/2)

template <int n, int d>
struct F0_wkspc {
  double sh;
  double u0;
  double du[d];
  double u_lam;
  double l;
  double q;
  double p_lam[n];
  double p0[n];
  double dP[d][n];
  double h; // TODO: delete?
  double theta; // TODO: delete?
};

template <int d, char... ps>
struct F0_bv_wkspc;

template <int d>
struct F0_bv_wkspc<d, char, char, char> {
  double sh;
  double u0;
  double du[d];
  double u_lam;
  double l;
  double y[d]; // TODO: more sensical name
  double h; // TODO: delete?
  double theta; // TODO: delete?
};

template <int n, int d>
struct F1_wkspc {
  double s_hat;
  double stheta;
  double s[d + 1];
  double ds[d];
  double q;
  double l;
  double u_lam;
  double u0;
  double du[d];
  double p_lam[n];
  double p0[n];
  double dP[d][n];
  double h;
  double theta;
};

template <int d, char... ps>
struct F1_bv_wkspc;

template <int d>
struct F1_bv_wkspc<d, char, char, char> {
  double s_hat;
  double sh;
  double s0;
  double ds[d];
  double l;
  double u_lam;
  double u0;
  double du[d];
  double y[d]; // TODO: more sensical name
  double h; // TODO: delete?
  double theta; // TODO: delete?
};

template <int n, int d>
struct fac_wkspc {
  double p_fac[n];
  double s_fac;
  double tau0;
  double dtau[d];
  double tau_lam;
  double dPt_dot_p_fac[d];
  double p_lam_minus_p_fac[n];
  double l_fac_lam;
  double l_fac_lam_sq;
};

template <int n, int d>
struct F0_fac: public cost_func<F0_fac<n, d>, n, d>
{
  F0_fac(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[d]) const;
  void hess_impl(double d2f[__sym_mat_size(d)]) const;
  void set_lambda_impl(double const lambda[d]);
  void set_args(double const u[d + 1], double s_hat,
                double const s[d + 1], double const p[d + 1][n],
                double const p_fac[n], double s_fac);
};  

template <int n, int d>
struct F1_fac: public cost_func<F1_fac<n, d>, n, d>
{
  F1_fac(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[d]) const;
  void hess_impl(double d2f[__sym_mat_size(d)]) const;
  void set_lambda_impl(double const lambda[d]);
  void set_args(double const u[d + 1], double s_hat,
                double const s[d + 1], double const p[d + 1][n],
                double const p_fac[n], double s_fac);
};

template <class derived, char p0, char p1, char p2, int d>
struct cost_func_bv {
  inline void eval(double & f) const {
    static_cast<derived const *>(this)->eval_impl(f);
  }

  inline void grad(double df[d]) const {
    static_cast<derived const *>(this)->grad_impl(df);
  }

  inline void hess(double d2f[__sym_mat_size(d)]) const {
    static_cast<derived const *>(this)->hess_impl(d2f);
  }

  void lag_mult(double const lambda[d], double * mu, int * k);

  inline void set_lambda(double const lambda[d]) {
    static_cast<derived *>(this)->set_lambda_impl(lambda);
  }

  inline void set_args(double const u[d + 1], double s_hat,
                       double const s[d + 1]) {
    static_cast<derived *>(this)->set_args_impl(u, s_hat, s);
  }
};

// template <char p0, char p1, char p2, int d>
// struct F0_bv: public cost_func_bv<F0_bv<p0, p1, p2, d>, p0, p1, p2, d> {
//   F0_bv(double h, double theta): _h {h}, _theta {theta} {}
//   void eval_impl(double & f) const;
//   void grad_impl(double df[d]) const;
//   void hess_impl(double d2f[__sym_mat_size(d)]) const;
//   void set_lambda_impl(double const lambda[d]);
//   void set_args_impl(double const u[d + 1], double s_hat,
//                      double const s[d + 1]);
// };

template <char p0, char p1, char p2>
struct F0_bv<p0, p1, p2, 2>:
  public cost_func_bv<F0_bv<p0, p1, p2, 2>, p0, p1, p2, 2>
{
  F0_bv(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[2]) const;
  void hess_impl(double d2f[3]) const;
  void set_lambda_impl(double const lambda[2]);
  void set_args_impl(double const u[3], double s_hat,
                     double const s[3]);
};

// template <char p0, char p1, char p2, int d>
// struct F1_bv: public cost_func_bv<F1_bv<p0, p1, p2, d>, p0, p1, p2, d> {
//   F1_bv(double h, double theta): _h {h}, _theta {theta} {}
//   void eval_impl(double & f) const;
//   void grad_impl(double df[d]) const;
//   void hess_impl(double d2f[__sym_mat_size(d)]) const;
//   void set_lambda_impl(double const lambda[d]);
//   void set_args_impl(double const u[d + 1], double s_hat,
//                      double const s[d + 1]);
// EIKONAL_PROTECTED:
// };

// Specialization for d = 2
template <char p0, char p1, char p2>
struct F1_bv<p0, p1, p2, 2>:
  public cost_func_bv<F1_bv<p0, p1, p2, 2>, p0, p1, p2, 2>
{
  F1_bv(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[2]) const;
  void hess_impl(double d2f[3]) const;
  void set_lambda_impl(double const lambda[2]);
  void set_args_impl(double const u[3], double s_hat,
                     double const s[3]);
EIKONAL_PROTECTED:
  double _s_hat;
  double _sh;
  double _s0;
  double _ds[2];
  double _l;
  double _u_lam;
  double _u0;
  double _du[2];
  double _y[2];
  double _h;
  double _theta;
};

#include "cost_funcs.impl.hpp"

#undef __sym_mat_size
#undef __check

#endif // __COST_FUNCS_HPP__
