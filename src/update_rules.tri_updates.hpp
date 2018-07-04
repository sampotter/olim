#ifndef __UPDATE_RULES_TRI_UPDATES_HPP__
#define __UPDATE_RULES_TRI_UPDATES_HPP__

#include <type_traits>

#include "common.macros.hpp"
#include "cost_funcs.hpp"
#include "update_rules.utils.hpp"

namespace update_rules {

struct mp0_tri_updates
{
  template <char p0, char p1>
  update_info<1> tri(
    double u0, double u1, double s, double s0, double s1, double h,
    ffvec<p0>, ffvec<p1>, double tol = EPS(double)) const;

  template <int dim>
  update_info<1> tri(
    double u0, double u1, double s, double s0, double s1, double h,
    double const * p0, double const * p1, double const * p_fac, double s_fac,
    double tol = EPS(double)) const;

  template <int d>
  update_info<1> tri(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h, double tol = EPS(double)) const;
};

struct mp1_tri_updates
{
  template <char p0, char p1>
  update_info<1> tri(
    double u0, double u1, double s, double s0, double s1, double h,
    ffvec<p0>, ffvec<p1>, double tol = EPS(double)) const;

  template <int dim>
  update_info<1> tri(
    double u0, double u1, double s, double s0, double s1, double h,
    double const * p0, double const * p1, double const * p_fac, double s_fac,
    double tol = EPS(double)) const;

  template <int d>
  update_info<1> tri(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h, double tol = EPS(double)) const;
};

struct rhr_tri_updates
{
  template <char p0, char p1>
  update_info<1> tri(
    double u0, double u1, double s, double s0, double s1, double h,
    ffvec<p0>, ffvec<p1>, double tol = EPS(double)) const;

  template <int dim>
  update_info<1> tri(
    double u0, double u1, double s, double s0, double s1, double h,
    double const * p0, double const * p1, double const * p_fac, double s_fac,
    double tol = EPS(double)) const;

  template <int d>
  update_info<1> tri(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h, double tol = EPS(double)) const;
};

}

#include "update_rules.tri_updates.impl.hpp"

#endif // __UPDATE_RULES_TRI_UPDATES_HPP__
