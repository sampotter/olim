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

/**
 * An explanation of the template arguments:
 * - derived: the effective derived subclass (cost_func uses CRTP)
 * - n: the dimension of the ambient space (where each p_i lives)
 * - d: the dimension of the domain of the cost function
 */
template <class derived, int n, int d>
struct cost_func {
  inline void eval(double & f) const {
    static_cast<derived const *>(this)->eval_impl(f);
    __check(f);
  }

  inline void grad(double df[d]) const {
    static_cast<derived const *>(this)->grad_impl(df);
    for (int i = 0; i < d; ++i) {
      __check(df[i]);
    }
  }

  inline void hess(double d2f[__sym_mat_size(d)]) const {
    static_cast<derived const *>(this)->hess_impl(d2f);
    for (int i = 0; i < __sym_mat_size(d); ++i) {
      __check(d2f[i]);
    }
  }

  void lag_mult(double const lambda[d], double * mu, int * k);

#if COLLECT_STATS
  inline bool degenerate_lambda() const {
    double lam0 = 0;
    for (int i = 0; i < d; ++i) {
      lam0 += _lam[i];
      if (_lam[i] < EPS(double)) return true;
    }
    return lam0 > 1 - EPS(double);
  }
#endif

  inline void set_lambda(double const lambda[d]) {
#if COLLECT_STATS
    for (int i = 0; i < d; ++i) _lam[i] = lambda[i];
#endif
    static_cast<derived *>(this)->set_lambda_impl(lambda);
  }

  inline void set_args(double const u[d + 1],
                       double s_hat,
                       double const s[d + 1],
                       double const p[d + 1][n]) {
    static_cast<derived *>(this)->set_args_impl(u, s_hat, s, p);
  }

#if COLLECT_STATS
EIKONAL_PROTECTED:
  double _lam[d];
#endif
};

template <int n, int d>
struct F0: public cost_func<F0<n, d>, n, d> {
  F0(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[d]) const;
  void hess_impl(double d2f[__sym_mat_size(d)]) const;
  void set_lambda_impl(double const lambda[d]);
  void set_args_impl(double const u[d + 1], double s_hat,
                     double const s[d + 1], double const p[d + 1][n]);
EIKONAL_PROTECTED:
  double _sh;
  double _u0;
  double _du[d];
  double _u_lam;
  double _l;
  double _q;
  double _p_lam[n];
  double _p0[n];
  double _dP[d][n];
  double _h;
  double _theta;
};

template <int n, int d>
struct F1: public cost_func<F1<n, d>, n, d> {
  F1(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[d]) const;
  void hess_impl(double d2f[__sym_mat_size(d)]) const;
  void set_lambda_impl(double const lambda[d]);
  void set_args_impl(double const u[d + 1], double s_hat,
                     double const s[d + 1], double const p[d + 1][n]);
EIKONAL_PROTECTED:
  double _s_hat;
  double _stheta;
  double _s[d + 1];
  double _ds[d];
  double _q;
  double _l;
  double _u_lam;
  double _u0;
  double _du[d];
  double _p_lam[n];
  double _p0[n];
  double _dP[d][n];
  double _h;
  double _theta;
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
EIKONAL_PROTECTED:
  double _sh;
  double _u0;
  double _du[d];
  double _u_lam;
  double _l;
  double _q;
  double _p_lam[n];
  double _p0[n];
  double _dP[d][n];
  double _h;
  double _theta;
  double _p_fac[n];
  double _s_fac;
  double _tau0;
  double _dtau[d];
  double _tau_lam;
  double _dPt_dot_p_fac[d];
  double _p_lam_minus_p_fac[n];
  double _l_fac_lam;
  double _l_fac_lam_sq;
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
EIKONAL_PROTECTED:
  double _s_hat;
  double _stheta;
  double _s[d + 1];
  double _ds[d];
  double _q;
  double _l;
  double _u_lam;
  double _u0;
  double _du[d];
  double _p_lam[n];
  double _p0[n];
  double _dP[d][n];
  double _h;
  double _theta;
  double _p_fac[n];
  double _s_fac;
  double _tau0;
  double _dtau[d];
  double _tau_lam;
  double _dPt_dot_p_fac[d];
  double _p_lam_minus_p_fac[n];
  double _l_fac_lam;
  double _l_fac_lam_sq;
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

#if COLLECT_STATS
  inline bool degenerate_lambda() const {
    double lam0 = 0;
    for (int i = 0; i < d; ++i) {
      lam0 += _lam[i];
      if (_lam[i] < EPS(double)) return true;
    }
    return lam0 > 1 - EPS(double);
  }
#endif

  inline void set_lambda(double const lambda[d]) {
#if COLLECT_STATS
    for (int i = 0; i < d; ++i) _lam[i] = lambda[i];
#endif
    static_cast<derived *>(this)->set_lambda_impl(lambda);
  }

  inline void set_args(double const u[d + 1], double s_hat,
                       double const s[d + 1]) {
    static_cast<derived *>(this)->set_args_impl(u, s_hat, s);
  }

#if COLLECT_STATS
EIKONAL_PROTECTED:
  double _lam[d];
#endif
};

template <char p0, char p1, char p2, int d>
struct F0_bv: public cost_func_bv<F0_bv<p0, p1, p2, d>, p0, p1, p2, d> {
  F0_bv(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[d]) const;
  void hess_impl(double d2f[__sym_mat_size(d)]) const;
  void set_lambda_impl(double const lambda[d]);
  void set_args_impl(double const u[d + 1], double s_hat,
                     double const s[d + 1]);
  // TODO: can probably delete this
EIKONAL_PROTECTED:
  double _sh;
  double _u0;
  double _du[d];
  double _u_lam;
  double _l;
  double _y[d];
  double _h;
  double _theta;
};

// Specialization for d = 2
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
EIKONAL_PROTECTED:
  double _sh;
  double _u0;
  double _du[2];
  double _u_lam;
  double _l;
  double _y[2];
  double _h;
  double _theta;
};

template <char p0, char p1, char p2, int d>
struct F1_bv: public cost_func_bv<F1_bv<p0, p1, p2, d>, p0, p1, p2, d> {
  F1_bv(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[d]) const;
  void hess_impl(double d2f[__sym_mat_size(d)]) const;
  void set_lambda_impl(double const lambda[d]);
  void set_args_impl(double const u[d + 1], double s_hat,
                     double const s[d + 1]);
EIKONAL_PROTECTED:
  double _s_hat;
  double _sh;
  double _s0;
  double _ds[d];
  double _l;
  double _u_lam;
  double _u0;
  double _du[d];
  double _y[d];
  double _h;
  double _theta;
};

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
