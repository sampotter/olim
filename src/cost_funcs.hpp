#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

template <int d>
struct F0 {
  F0(double h, double theta): _h {h}, _theta {theta} {}
  void eval(double & f) const;
  void grad(double df[d]) const;
  void hess(double d2f[d][d]) const;
  void set_lambda(double const lambda[d]);
  void set_args(double const u[d + 1], double s_hat, double const s[d + 1],
                double const p[d + 1][d + 1]);
private:
  double _sh;
  double _u0;
  double _du[d];
  double _u_lam;
  double _l;
  double _q;
  double _p_lam[d + 1];
  double _p0[d + 1];
  double _dP[d][d + 1];
  double _h;
  double _theta;
};

template <int d>
struct F1 {
  F1(double h, double theta): _h {h}, _theta {theta} {}
  void eval(double & f) const;
  void grad(double df[d]) const;
  void hess(double d2f[d][d]) const;
  void set_lambda(double const lambda[d]);
  void set_args(double const u[d + 1], double s_hat, double const s[d + 1],
                double const p[d + 1][d + 1]);
private:
  double _s_hat;
  double _stheta;
  double _s[d + 1];
  double _ds[d];
  double _q;
  double _l;
  double _u_lam;
  double _u0;
  double _du[d];
  double _p_lam[d + 1];
  double _p0[d + 1];
  double _dP[d][d + 1];
  double _h;
  double _theta;
};

#endif // __COST_FUNCS_HPP__
