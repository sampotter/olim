#ifndef __UPDATE_RULES_TETRA_UPDATES_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_HPP__

#include "cost_funcs.hpp"
#include "update_rules.utils.hpp"

namespace update_rules {

template <class derived>
struct tetra_updates
{
  template <char p0, char p1, char p2>
  update_info<2> tetra(
    double u0, double u1, double u2, double s,
    double s0, double s1, double s2, double h,
    ffvec<p0>, ffvec<p1>, ffvec<p2>) const;

  update_info<2> tetra(
    double u0, double u1, double u2, double s,
    double s0, double s1, double s2, double h,
    double const * p0, double const * p1, double const * p2,
    double const * p_fac, double s_fac) const;

  update_info<2> tetra(
    double const * p0, double const * p1, double const * p2,
    double u0, double u1, double u2, double s,
    double s0, double s1, double s2, double h) const;
};

struct mp0_tetra_updates: tetra_updates<mp0_tetra_updates> {
  template <int n, int d>
  using cost_func = F0<n, d>;

  template <char p0, char p1, char p2>
  using cost_func_bv = F0_bv<p0, p1, p2, 2>;

  template <int n, int d>
  using factored_cost_func = F0_fac<n, d>;

  inline double theta() const { return 0.5; }
};

struct mp1_tetra_updates: tetra_updates<mp1_tetra_updates> {
  template <int n, int d>
  using cost_func = F1<n, d>;

  template <char p0, char p1, char p2>
  using cost_func_bv = F1_bv<p0, p1, p2, 2>;

  template <int n, int d>
  using factored_cost_func = F1_fac<n, d>;

  inline double theta() const { return 0.5; }
};

struct rhr_tetra_updates: tetra_updates<rhr_tetra_updates> {
  template <int n, int d>
  using cost_func = F0<n, d>;

  template <char p0, char p1, char p2>
  using cost_func_bv = F0_bv<p0, p1, p2, 2>;

  template <int n, int d>
  using factored_cost_func = F0_fac<n, d>;

  inline double theta() const { return 0.0; }
};

}

#include "update_rules.tetra_updates.impl.hpp"

#endif // __UPDATE_RULES_TETRA_UPDATES_HPP__
