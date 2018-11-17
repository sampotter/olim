#ifndef __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#if PRINT_UPDATES
#    include <cstdio>
#endif
#include <type_traits>

#include "common.defs.hpp"
#include "common.macros.hpp"
#include "cost_funcs.hpp"
#include "numopt.hpp"
#include "update_rules.tetra_updates.util.hpp"

template <cost_func F>
update_info<2> update_rules::tetra<F>::operator()(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h) const
{
  struct {
    inline void eval(double & f) const { ::eval(w, f); }
    inline void grad(double * df) const { ::grad(w, df); }
    inline void hess(double * d2f) const { ::hess(w, d2f); }
    inline void set_lambda(double * lam) {
      ::set_lambda<F, 3>(w, p0, p1, p2, lam);
    }
    F_wkspc<F, 2> w;
    double const * p0, * p1, * p2;
  } func;

  func.p0 = p0;
  func.p1 = p1;
  func.p2 = p2;

  set_args<F, 3>(func.w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);

  update_info<2> update;
  bool error;
  sqp_bary<decltype(func), 3, 2> sqp;
  sqp(func, update.lambda, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    // TODO: check and see if this can be optimized at all to avoid
    // redundant calculations with set_args call above
    F_wkspc<MP1, 2> w;
    set_args<MP1, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h);
    set_lambda<MP1, 3>(w, p0, p1, p2, update.lambda);
    eval(w, update.value);
  } else {
    // TODO: we're doing an unnecessary eval here: we could reorganize
    // things so that we're using the most recent eval done by
    // sqp... i.e., have sqp write over update.value internally
    eval(func.w, update.value);
  }

#if PRINT_UPDATES
  printf("tetra(u0 = %g, u1 = %g, u2 = %g, s = %g, "
         "s0 = %g, s1 = %g, s2 = %g, h = %g) -> %g\n",
         u0, u1, u2, s, s0, s1, s2, h, update.value);
#endif

  return update;
}

template <cost_func F>
update_info<2> update_rules::tetra<F>::operator()(
  double const * p0, double const * p1, double const * p2,
  double u0, double u1, double u2, double s,
  double s0, double s1, double s2, double h,
  double const * p_fac, double s_fac) const  
{
  struct {
    inline void eval(double & f) const { ::eval(w, f); }
    inline void grad(double * df) const { ::grad(w, df); }
    inline void hess(double * d2f) const { ::hess(w, d2f); }
    inline void set_lambda(double const * lam) {
      ::set_lambda<3>(w, p0, p1, p2, p_fac, lam);
    }
    F_fac_wkspc<F, 2> w;
    double const * p0, * p1, * p2, * p_fac;
  } func;

  func.p0 = p0;
  func.p1 = p1;
  func.p2 = p2;
  func.p_fac = p_fac;
  
  set_args<F, 3>(func.w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, p_fac, s_fac);
  
  update_info<2> update;
  bool error;
  sqp_bary<decltype(func), 3, 2> sqp;
  sqp(func, update.lambda, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    // TODO: check and see if this can be optimized at all to avoid
    // redundant calculations with set_args call above
    F_fac_wkspc<MP1, 2> w;
    set_args<MP1, 3>(w, p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, p_fac, s_fac);
    set_lambda<3>(w, p0, p1, p2, p_fac, update.lambda);
    eval(w, update.value);
  } else {
    eval(func.w, update.value);
  }

#if PRINT_UPDATES
#  error Not implemented yet
#endif

  return update;
}

template <cost_func F, char p0, char p1, char p2>
update_info<2>
update_rules::tetra_bv<F, p0, p1, p2>::operator()(F_wkspc<F, 2> const & w) const
{
  update_info<2> update;
  bool error;
  sqp_bary<F_wkspc<F, 2>, 3, 2> sqp;
  sqp(w, update.lambda, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    assert(false); // fix later
    // TODO: check and see if this can be optimized at all to avoid
    // redundant calculations with set_args call above
    // F_wkspc<MP1, 2> w;
    // set_args<MP1, 3, p0, p1, p2>(w, u0, u1, u2, s, s0, s1, s2, h);
    // set_lambda<MP1, 3, p0, p1, p2>(w, update.lambda);
    // eval(w, update.value);
  } else {
    // TODO: I think we should actually be okay to remove this---try
    // doing so once tests are stabilized after this big change
    eval(w, update.value);
  }

#if PRINT_UPDATES
  printf("tetra<%d, %d, %d>(u0 = %g, u1 = %g, u2 = %g, s = %g, "
         "s0 = %g, s1 = %g, s2 = %g, h = %g) -> %g\n",
         p0, p1, p2, u0, u1, u2, s, s0, s1, s2, h, update.value);
#endif

  return update;
}

template <cost_func F, char p0, char p1, char p2>
update_info<2>
update_rules::tetra_bv<F, p0, p1, p2>::operator()(F_fac_wkspc<F, 2> const & fw) const
{
  // struct {
  //   inline void eval(double & f) const { eval(w, fw, f); }
  //   inline void grad(double * df) const { grad(w, fw, df); }
  //   inline void hess(double * d2f) const { hess(w, fw, d2f); }
  //   inline void set_lambda(double const * lam) const {
  //     ::set_lambda<3, p0, p1, p2>(w, lam);
  //     ::set_lambda<3, p0, p1, p2>(fw, lam);
  //   }
  //   F_wkspc<F, 2> w;
  //   fac_wkspc<2> fw;
  // } func;

  // set_args<F, 3, p0, p1, p2>(func.w, u0, u1, u2, s, s0, s1, s2, h);
  // set_args<F, 3, p0, p1, p2>(func.fw, u0, u1, u2, s, s0, s1, s2, h, p_fac, s_fac);
  
  update_info<2> update;
  bool error;
  sqp_bary<F_fac_wkspc<F, 2>, 3, 2> sqp;
  sqp(fw, update.lambda, &error);
  assert(!error);

  if (F == cost_func::mp0) {
    assert(false); // fix later
    // TODO: check and see if this can be optimized at all to avoid
    // redundant calculations with set_args call above
    // F_wkspc<MP1, 2> w;
    // set_args<MP1, 3, p0, p1, p2>(w, u0, u1, u2, s, s0, s1, s2, h);
    // set_lambda<MP1, 3, p0, p1, p2>(w, update.lambda);
    // eval(w, func.fw, update.value);
  } else {
    // TODO: can probably remove
    eval(fw, update.value);
  }

#if PRINT_UPDATES
#  error Not implemented yet
#endif

  return update;
}

#undef __theta

#endif // __UPDATE_RULES_TETRA_UPDATES_IMPL_HPP__
