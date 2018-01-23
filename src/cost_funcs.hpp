#ifndef __COST_FUNCS_HPP__
#define __COST_FUNCS_HPP__

// TODO: this is an awful hack that we're going to use for now---we
// obviously don't want to use virtual functions in a hot path like
// this. We're just going this route until we take care of some other
// things that are higher priority.

/**
 * An explanation of the template arguments:
 * - derived: the effective derived subclass (this implementation uses
 *   CRTP)
 * - n: the dimension of the ambient space
 * - d: the affine dimension of the simplex (i.e. domain of the cost
 *   function)
 */
template <class derived, int n, int d>
struct cost_func {
  inline void eval(double & f) const {
    static_cast<derived const *>(this)->eval_impl(f);
  }

  inline void grad(double df[d]) const {
    static_cast<derived const *>(this)->grad_impl(df);
  }

  inline void hess(double d2f[d][d]) const {
    static_cast<derived const *>(this)->hess_impl(d2f);
  }

  inline void set_lambda(double const lambda[d]) {
    static_cast<derived *>(this)->set_lambda_impl(lambda);
  }

  inline void set_args(double const u[d + 1],
                       double s_hat,
                       double const s[d + 1],
                       double const p[d + 1][n]) {
    static_cast<derived *>(this)->set_args_impl(u, s_hat, s, p);
  }
};

template <int n, int d>
struct F0: public cost_func<F0<n, d>, n, d> {
  F0(double h, double theta): _h {h}, _theta {theta} {}
  void eval_impl(double & f) const;
  void grad_impl(double df[d]) const;
  void hess_impl(double d2f[d][d]) const;
  void set_lambda_impl(double const lambda[d]);
  void set_args_impl(double const u[d + 1], double s_hat,
                     double const s[d + 1], double const p[d + 1][n]);
private:
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
  void hess_impl(double d2f[d][d]) const;
  void set_lambda_impl(double const lambda[d]);
  void set_args_impl(double const u[d + 1], double s_hat,
                     double const s[d + 1], double const p[d + 1][n]);
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
  double _p_lam[n];
  double _p0[n];
  double _dP[d][n];
  double _h;
  double _theta;
};

#endif // __COST_FUNCS_HPP__
