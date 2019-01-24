#ifndef __UPDATES_TRI_HPP__
#define __UPDATES_TRI_HPP__

// TODO: remove me
#include <type_traits>

#include "common.hpp"
#include "cost_funcs.hpp"
#include "updates.common.hpp"
#include "vec.hpp"

// TODO: this is a bit of a mess right now---we would like to collapse
// the implementations like we've done for line and tetra, but that
// will have to wait until later.

namespace updates {

template <cost_func F, int n> struct tri {};

template <int n>
struct tri<MP0, n>
{
  info<1> operator()(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h) const;
    
  info<1> operator()(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h,
    double const * p_fac, double s_fac) const;
};

template <int n>
struct tri<MP1, n>
{
  info<1> operator()(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h) const;
    
  info<1> operator()(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h,
    double const * p_fac, double s_fac) const;
};

template <int n>
struct tri<RHR, n>
{
  info<1> operator()(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h) const;
    
  info<1> operator()(
    double const * p0, double const * p1, double u0, double u1,
    double s, double s0, double s1, double h,
    double const * p_fac, double s_fac) const;
};

template <cost_func F, int n, int p0, int p1> struct tri_bv {};

template <int n, int p0, int p1>
struct tri_bv<MP0, n, p0, p1>
{
  info<1> operator()(
    double u0, double u1, double s, double s0, double s1, double h) const;
};

template <int n, int p0, int p1>
struct tri_bv<MP1, n, p0, p1>
{
  info<1> operator()(
    double u0, double u1, double s, double s0, double s1, double h) const;
};

template <int n, int p0, int p1>
struct tri_bv<RHR, n, p0, p1>
{
  info<1> operator()(
    double u0, double u1, double s, double s0, double s1, double h) const;
};

}

#include "updates.tri.impl.hpp"

#endif // __UPDATES_TRI_HPP__
