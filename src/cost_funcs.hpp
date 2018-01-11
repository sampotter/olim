#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

// TODO: this is an awful hack that we're going to use for now---we
// obviously don't want to use virtual functions in a hot path like
// this. We're just going this route until we take care of some other
// things that are higher priority.

template <int d>
struct cost_func {
  virtual void eval(double & f) const = 0;
  virtual void grad(double df[d]) const = 0;
  virtual void hess(double d2f[d][d]) const = 0;
  virtual void set_lambda(double const lambda[d]) = 0;
  virtual void set_args(double const u[d + 1], double s_hat,
                        double const s[d + 1],
                        double const p[d + 1][d + 1]) = 0;
};

template <int d>
struct F0: public cost_func<d> {
  F0(double h, double theta): _h {h}, _theta {theta} {}
  virtual void eval(double & f) const;
  virtual void grad(double df[d]) const;
  virtual void hess(double d2f[d][d]) const;
  virtual void set_lambda(double const lambda[d]);
  virtual void set_args(double const u[d + 1], double s_hat,
                        double const s[d + 1], double const p[d + 1][d + 1]);
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
struct F1: public cost_func<d> {
  F1(double h, double theta): _h {h}, _theta {theta} {}
  virtual void eval(double & f) const;
  virtual void grad(double df[d]) const;
  virtual void hess(double d2f[d][d]) const;
  virtual void set_lambda(double const lambda[d]);
  virtual void set_args(double const u[d + 1], double s_hat,
                        double const s[d + 1], double const p[d + 1][d + 1]);
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
