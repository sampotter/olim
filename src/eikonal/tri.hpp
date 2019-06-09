#pragma once

// TODO: remove me
#include <type_traits>

#include "../type.h"

#include "../common.hpp"
#include "../update_info.hpp"
#include "../vec.hpp"

// TODO: this is a bit of a mess right now---we would like to collapse
// the implementations like we've done for line and tetra, but that
// will have to wait until later.

namespace updates {

template <cost_func F, int n> struct tri {};

template <int n>
struct tri<MP0, n>
{
  update_info<1> operator()(
    vec<double, n> const & p0, vec<double, n> const & p1, double u0, double u1,
    double s, double s0, double s1, double h) const;

  update_info<1> operator()(
    vec<double, n> const & p0, vec<double, n> const & p1, double u0, double u1,
    double s, double s0, double s1, double h,
    vec<double, n> const & p_fac, double s_fac) const;
};

template <int n>
struct tri<MP1, n>
{
  update_info<1> operator()(
    vec<double, n> const & p0, vec<double, n> const & p1, double u0, double u1,
    double s, double s0, double s1, double h) const;

  update_info<1> operator()(
    vec<double, n> const & p0, vec<double, n> const & p1, double u0, double u1,
    double s, double s0, double s1, double h,
    vec<double, n> const & p_fac, double s_fac) const;
};

template <int n>
struct tri<RHR, n>
{
  update_info<1> operator()(
    vec<double, n> const & p0, vec<double, n> const & p1, double u0, double u1,
    double s, double s0, double s1, double h) const;

  update_info<1> operator()(
    vec<double, n> const & p0, vec<double, n> const & p1, double u0, double u1,
    double s, double s0, double s1, double h,
    vec<double, n> const & p_fac, double s_fac) const;
};

template <cost_func F, int n, int p0, int p1> struct tri_bv {};

template <int n, int p0, int p1>
struct tri_bv<MP0, n, p0, p1>
{
  update_info<1> operator()(
    double u0, double u1, double s, double s0, double s1, double h) const;
};

template <int n, int p0, int p1>
struct tri_bv<MP1, n, p0, p1>
{
  update_info<1> operator()(
    double u0, double u1, double s, double s0, double s1, double h) const;
};

template <int n, int p0, int p1>
struct tri_bv<RHR, n, p0, p1>
{
  update_info<1> operator()(
    double u0, double u1, double s, double s0, double s1, double h) const;
};

}

#include "tri.impl.hpp"
