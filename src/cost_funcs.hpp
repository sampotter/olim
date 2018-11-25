#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

#include <cassert>
#include <cmath>
#include <type_traits>

#include "common.macros.hpp"

enum class cost_func {mp0, mp1, rhr};

constexpr auto MP0 = cost_func::mp0;
constexpr auto MP1 = cost_func::mp1;
constexpr auto RHR = cost_func::rhr;

// template <cost_func F>
// constexpr double theta = F == MP0 || F == MP1 ? 0.5 : 0;

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

template <cost_func F, int d>
using F_fac_wkspc = std::conditional_t<
  F == MP1, F1_wkspc<d, fac_wkspc<d>>, F0_wkspc<d, fac_wkspc<d>>>;

// TODO: below are specializations for d = 2

template <cost_func F, int n, char p0, char p1, char p2>
void set_args(F_wkspc<F, 2> & w, double u0, double u1, double u2,
              double s, double s0, double s1, double s2, double h);

template <cost_func F, int n>
void set_args(F_wkspc<F, 2> & w,
              double const * p0, double const * p1, double const * p2,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h);

template <cost_func F, int n>
void set_args(F_fac_wkspc<F, 2> & w,
              double const * p0, double const * p1, double const * p2,
              double u0, double u1, double u2, double s,
              double s0, double s1, double s2, double h,
              double const * p_fac, double s_fac);

template <cost_func F, int n, char p0, char p1, char p2>
void set_lambda(F_wkspc<F, 2> & w, double const * lam);

template <cost_func F, int n>
void set_lambda(F_wkspc<F, 2> & w, double const * p0, double const * p1,
                double const * p2, double const * lam);

template <int n>
void set_lambda(F0_wkspc<2, fac_wkspc<2>> & w,
                double const * p0, double const * p1, double const * p2,
                double const * p_fac, double const * lam);

template <int n>
void set_lambda(F1_wkspc<2, fac_wkspc<2>> & w,
                double const * p0, double const * p1, double const * p2,
                double const * p_fac, double const * lam);

// template <int d>
// void eval(eval_wkspc<d> const & w, double & f);

// template <int d>
// void eval(fac_wkspc<d> const & w, double & f);

template <int d>
void grad(F0_wkspc<d> const & w, double * df);

template <int d>
void grad(F0_wkspc<d, fac_wkspc<d>> const & w, double * df);

template <int d>
void grad(F1_wkspc<d> const & w, double * df);

template <int d>
void grad(F1_wkspc<d, fac_wkspc<d>> const & w, double * df);

template <int d>
void hess(F0_wkspc<d> const & w, double * d2f);

template <int d>
void hess(F0_wkspc<d, fac_wkspc<d>> const & w, double * d2f);

template <int d>
void hess(F1_wkspc<d> const & w, double * d2f);

template <int d>
void hess(F1_wkspc<d, fac_wkspc<d>> const & w, double * d2f);

#include "cost_funcs.impl.hpp"

template <cost_func F, char... ps>
struct cost_functor_bv;

template <cost_func F, int n, char p0, char p1, char p2>
struct cost_functor_bv<F, n, p0, p1, p2>
{
  using wkspc = F_wkspc<F, 2>;
  cost_functor_bv(F_wkspc<F, 2> & w): w {w} {}
  void set_lambda(double const * lam) { ::set_lambda<F, n, p0, p1, p2>(w, lam); }
  void eval(double & f) { ::eval(w, f); }
  void grad(double * df) { ::grad(w, df); }
  void hess(double * d2f) { ::hess(w, d2f); }
  wkspc & w;
};

#endif // __COST_FUNCS_HPP__
