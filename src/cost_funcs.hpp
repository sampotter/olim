#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

constexpr int sym_mat_size(int d) {
  return ((d + 1)*d)/2;
}

template <int d>
void lagmults(double const * lam, double const * df, double const * d2f,
              double * mu, int * k);

struct eval_wkspc {
  double sh;
  double u_lam;
  double l_lam;
};

template <int d>
struct F0_wkspc: public eval_wkspc {
  double du[d];
  double dPt_nu_lam[d];
  double dPt_dP[sym_mat_size(d)];
};

// template <int n, int d>
// struct F0_wkspc: public eval_wkspc {
  // double u0;
  // double du[d];
  // double q;
  // double p_lam[n];
  // double p0[n];
  // double dP[d][n];
// };

// template <int d, class... ps>
// struct F0_bv_wkspc;

// template <int d>
// struct F0_bv_wkspc<d, char, char, char>: public eval_wkspc {
  // double u0;
  // double du[d];
  // double y[d]; // TODO: more sensical name
// };

template <int d>
struct F1_wkspc: public eval_wkspc {
  double du[d];
  double dPt_nu_lam[d];
  double dPt_dP[sym_mat_size(d)];
  double theta_h_ds[d];
};

// template <int n, int d>
// struct F1_wkspc: public eval_wkspc {
  // double s_hat;
  // double stheta;
  // double s[d + 1];
  // double ds[d];
  // double q;
  // double u0;
  // double du[d];
  // double p_lam[n];
  // double p0[n];
  // double dP[d][n];
// };

// template <int d, class... ps>
// struct F1_bv_wkspc;

// template <int d>
// struct F1_bv_wkspc<d, char, char, char>: public eval_wkspc {
  // double s_hat;
  // double s0;
  // double ds[d];
  // double u0;
  // double du[d];
  // double y[d]; // TODO: more sensical name
// };

struct fac_eval_wkspc {
  double tau_lam;
  double sh_fac;
  double l_fac_lam;
};

template <int d>
struct fac_wkspc: public fac_eval_wkspc {
  double dPt_nu_fac_lam[d];
};

// template <int n, int d>
// struct fac_wkspc: public fac_eval_wkspc {
  // double p_fac[n];
  // double tau0;
  // double dtau[d];
  // double dPt_dot_p_fac[d];
  // double p_lam_minus_p_fac[n];
  // double l_fac_lam_sq;
// };

void eval(eval_wkspc const & w, double & f);
void eval(eval_wkspc const & w, fac_eval_wkspc const & fw, double & f);

template <int d>
void grad(F0_wkspc<d> const & w, double * df);

template <int d>
void grad(F0_wkspc<d> const & w, fac_wkspc<d> const & fw, double * df);

template <int d>
void grad(F1_wkspc<d> const & w, double * df);

template <int d>
void grad(F1_wkspc<d> const & w, fac_wkspc<d> const & fw, double * df);

template <int d>
void hess(F0_wkspc<d> const & w, double * d2f);

template <int d>
void hess(F0_wkspc<d> const & w, fac_wkspc<d> const & fw, double * d2f);

template <int d>
void hess(F1_wkspc<d> const & w, double * d2f);

template <int d>
void hess(F1_wkspc<d> const & w, fac_wkspc<d> const & fw, double * d2f);

#include "cost_funcs.impl.hpp"

#endif // __COST_FUNCS_HPP__
