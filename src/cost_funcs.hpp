#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

template <int d>
void lagmults(double const * lam, double const * df, double const * d2f,
              double * mu, int * k);

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

template <int d, class... ps>
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

template <int d, class... ps>
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

template <class wkspc>
void eval(wkspc const & w, double & f);

template <class wkspc>
void eval(wkspc const & w, fac_wkspc const & fw, double & f);

template <class wkspc>
void grad(wkspc const & w, double * df);

template <class wkspc>
void grad(wkspc const & w, fac_wkspc const & fw, double * df);

template <class wkspc>
void hess(wkspc const & w, double * d2f);

template <class wkspc>
void hess(wkspc const & w, fac_wkspc const & fw, double * d2f);

#include "cost_funcs.impl.hpp"

#endif // __COST_FUNCS_HPP__
