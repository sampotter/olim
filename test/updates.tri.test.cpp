#include <gtest/gtest.h>

#include <src/config.hpp>

#include "common.hpp"
#include "updates.tri.hpp"

using namespace updates;

#define P01 1
#define P10 2
#define P11 3

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

template <cost_func F>
info<1> tri11(double u0, double u1, double s, double s0, double s1, double h) {
  return tri_bv<F, 2, P01, P10>()(u0, u1, s, s0, s1, h);
}

template <cost_func F>
info<1> tri12(double u0, double u1, double s, double s0, double s1, double h) {
  return tri_bv<F, 2, P01, P11>()(u0, u1, s, s0, s1, h);
}

template <cost_func F>
info<1> tri13(double u0, double u1, double s, double s0, double s1, double h) {
  return tri_bv<F, 3, P001, P111>()(u0, u1, s, s0, s1, h);
}

template <cost_func F>
info<1> tri22(double u0, double u1, double s, double s0, double s1, double h) {
  return tri_bv<F, 3, P011, P101>()(u0, u1, s, s0, s1, h);
}

template <cost_func F>
info<1> tri23(double u0, double u1, double s, double s0, double s1, double h) {
  return tri_bv<F, 3, P011, P111>()(u0, u1, s, s0, s1, h);
}

/**
 * rhr tri11 tests
 */

TEST (updates_tri, rhr_tri11_basic_test) {
  double u = tri11<RHR>(0.0, 0.0, 1.0, 1.0, 1.0, 1.0).value;
  ASSERT_DOUBLE_EQ(u, sqrt2/2);
  u = tri11<RHR>(0.0, 1.0, 1.0, 1.0, 1.0, 1.0).value;
  ASSERT_DOUBLE_EQ(u, 1.0);
  u = tri11<RHR>(1.0, 0.0, 1.0, 1.0, 1.0, 1.0).value;
  ASSERT_DOUBLE_EQ(u, 1.0);
}

template <cost_func F>
void tri11_is_symmetric_with_constant_slowness() {
  double U0 = 0.0, U1 = 0.1, s = 1, s0 = 1, s1 = 1, h = 1, u, u1, u2;
  u = 0.7553367989832942;
  u1 = tri11<F>(U0, U1, s, s0, s1, h).value;
  u2 = tri11<F>(U1, U0, s, s1, s0, h).value;
  ASSERT_DOUBLE_EQ(u1, u);
  ASSERT_DOUBLE_EQ(u2, u);
  ASSERT_DOUBLE_EQ(u1, u2);
}

TEST (updates_tri, tri11_is_symmetric_with_constant_slowness) {
  tri11_is_symmetric_with_constant_slowness<MP0>();
  tri11_is_symmetric_with_constant_slowness<MP1>();
  tri11_is_symmetric_with_constant_slowness<RHR>();
}

template <cost_func F>
void tri11_is_symmetric_with_nonconstant_slowness(
  double U, double u0, double u1, double s, double s0, double s1, double h)
{
  double U1 = tri11<F>(u0, u1, s, s0, s1, h).value;
  double U2 = tri11<F>(u1, u0, s, s1, s0, h).value;
  ASSERT_DOUBLE_EQ(U1, U);
  ASSERT_DOUBLE_EQ(U2, U);
  ASSERT_DOUBLE_EQ(U1, U2);
}

TEST (updates_tri, tri11_is_symmetric_with_nonconstant_slowness) {
  double u0 = 0, u1 = 0.1, s = 1, s0 = 1.2, s1 = 1.1, h = 1;
  double u_mp0 = 0.809661416377184, u_mp1 = 0.809452960468600,
    u_rhr = 0.7553367989832942;
  tri11_is_symmetric_with_nonconstant_slowness<MP0>(u_mp0, u0, u1, s, s0, s1, h);
  tri11_is_symmetric_with_nonconstant_slowness<MP1>(u_mp1, u0, u1, s, s0, s1, h);
  tri11_is_symmetric_with_nonconstant_slowness<RHR>(u_rhr, u0, u1, s, s0, s1, h);
}

TEST (updates_tri, rhr_tri12_works) {
  double t, u0, u1, s = 1, s0 = 1, s1 = 1, h = 0.25;

  t = tri12<RHR>(1.0, 0.0, 1.0, 1.0, 1.0, 1.0).value;
  ASSERT_DOUBLE_EQ(t, sqrt2);

  ASSERT_DOUBLE_EQ(tri12<RHR>(0.0, 1.0, 1.0, 1.0, 1.0, 1.0).value, 1.0);

  u0 = 0.923028, u1 = 0.707107;
  ASSERT_NEAR(tri12<RHR>(u0, u1, s, s0, s1, h).value, 1.060660390593274, 1e-14);

  u0 = 0.852848, u1 = 0.707107;
  ASSERT_NEAR(tri12<RHR>(u0, u1, s, s0, s1, h).value, 1.05597249, 6.2e-9);

  u0 = 0.923028, u1 = 0.817579;
  ASSERT_NEAR(tri12<RHR>(u0, u1, s, s0, s1, h).value, 1.14970069, 1.6e-11);

  u0 = 0.883639, u1 = 0.817579;
  ASSERT_NEAR(tri12<RHR>(u0, u1, s, s0, s1, h).value, 1.12475323, 9.32e-9);

  u0 = 0.852848, u1 = 0.65974;
  ASSERT_NEAR(tri12<RHR>(u0, u1, s, s0, s1, h).value, 1.013293390593274, 1e-14);

  u0 = 0.883639, u1 = 0.65974;
  ASSERT_NEAR(tri12<RHR>(u0, u1, s, s0, s1, h).value, 1.013293390593274, 1e-14);
}

TEST (updates_tri, rhr_tri13_works) {
  ASSERT_DOUBLE_EQ(tri13<RHR>(1.0, 0.0, 1.0, 1.0, 1.0, 1.0).value, 1.707106781186547);
  ASSERT_DOUBLE_EQ(tri13<RHR>(0.0, 1.0, 1.0, 1.0, 1.0, 1.0).value, 1.0);
}

TEST (updates_tri, rhr_tri22_works) {
  ASSERT_DOUBLE_EQ(tri22<RHR>(0.0, 0.0, 1.0, 1.0, 1.0, 1.0).value, 1.224744871391589);
  ASSERT_DOUBLE_EQ(tri22<RHR>(0.0, 1.0, 1.0, 1.0, 1.0, 1.0).value, 1.414213562373095);
  ASSERT_DOUBLE_EQ(tri22<RHR>(1.0, 0.0, 1.0, 1.0, 1.0, 1.0).value, 1.414213562373095);

  double u0, u1, s = 1, s0 = 1, s1 = 1, h = 0.25;

  u0 = 0.65974, u1 = 0.817579;
  ASSERT_DOUBLE_EQ(tri22<RHR>(u0, u1, s, s0, s1, h).value, 1.012639677310786);

  u0 = 0.65974, u1 = 0.707107;
  ASSERT_DOUBLE_EQ(tri22<RHR>(u0, u1, s, s0, s1, h).value, 0.986849397845339);

  u0 = 0.707107, u1 = 0.817579;
  ASSERT_DOUBLE_EQ(tri22<RHR>(u0, u1, s, s0, s1, h).value, 1.053198553345643);
}

template <cost_func F>
void tri22_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, s0 = 1.0, s1 = 1.0, h = 1;
  double val01 = tri22<F>(u0, u1, s, s0, s1, h).value;
  double val10 = tri22<F>(u1, u0, s, s1, s0, h).value;
  ASSERT_DOUBLE_EQ(val01, val10);
}

TEST (updates_tri, tri22_is_symmetric_with_constant_slowness) {
  tri22_is_symmetric<MP0>();
  tri22_is_symmetric<MP1>();
  tri22_is_symmetric<RHR>();
}

TEST (updates_tri, rhr_tri23_works) {
  ASSERT_DOUBLE_EQ(tri23<RHR>(1.0, 0.0, 1.0, 1.0, 1.0, 1.0).value, 1.732050807568877);
  ASSERT_DOUBLE_EQ(tri23<RHR>(0.0, 1.0, 1.0, 1.0, 1.0, 1.0).value, 1.414213562373095);
}

TEST (updates_tri, mp1_tri11_works) {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.1, u1 = 0, s = 1, s0 = 1.2, s1 = 1.1, h = 0.9, T = 0.7309362364283433;
  ASSERT_NEAR(tri11<MP1>(u0, u1, s, s0, s1, h).value, T, 1e-13);
  ASSERT_NEAR(tri11<MP1>(u1, u0, s, s1, s0, h).value, T, 1e-13);

  u0 = 0, u1 = 0.1, s = 1, s0 = 1.1, s1 = 1.3, h = 1.2, T = 0.978772426584215;
  ASSERT_NEAR(tri11<MP1>(u0, u1, s, s0, s1, h).value, T, 1e-13);
  ASSERT_NEAR(tri11<MP1>(u1, u0, s, s1, s0, h).value, T, 1e-13);

  u0 = 1.0, u1 = 0.853553, s = 1.0, s0 = 1.0, s1 = 1.0, h = 0.5,
    T = 1.27266422607271;
  ASSERT_NEAR(tri11<MP1>(u0, u1, s, s0, s1, h).value, T, 1e-13);
  ASSERT_NEAR(tri11<MP1>(u1, u0, s, s1, s0, h).value, T, 1e-13);
}

TEST (updates_tri, mp1_tri12_works) {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.46487865669951, u1 = 0.4, s = 1, s0 = 1, s1 = 1, h = 0.2;
  T = 0.6540631166978595;
  ASSERT_DOUBLE_EQ(tri12<MP1>(u0, u1, s, s0, s1, h).value, T);
}

TEST (updates_tri, mp1_tri13_works) {
  double u0, u1, s, s0, s1, h;

  u0 = 2, u1 = 1, s = 1, s0 = 1, s1 = 1, h = 1;
  ASSERT_DOUBLE_EQ(tri13<MP1>(u0, u1, s, s0, s1, h).value, 2 + 1./sqrt(2));
}

TEST (updates_tri, mp1_tri22_works) {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.65974, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  T = 0.986849397845339;
  ASSERT_DOUBLE_EQ(tri22<MP1>(u0, u1, s, s0, s1, h).value, T);

  u0 = 0.59457611725233628;
  u1 = 0.56102682919361013;
  s = 0.0025050133959455545;
  s0 = 0.0025050133959455545;
  s1 = 0.2382400185837108;
  h = 0.5;
  T = 0.5963474292115725;
  ASSERT_DOUBLE_EQ(tri22<MP1>(u0, u1, s, s0, s1, h).value, T);
}

TEST (updates_tri, rhr_basic_factoring_test) {
  double u0, u1, s, s0, s1, h, s_fac;
  vec<double, 2> p0, p1, p_fac;

  u0 = u1 = s = s0 = s1 = h = s_fac = 1;

  p0[0] = 1;
  p0[1] = 0;

  p1[0] = 0;
  p1[1] = 1;

  p_fac[0] = 1;
  p_fac[1] = 1;

  auto update = tri<RHR, 2>()(p0, p1, u0, u1, s, s0, s1, h, p_fac, s_fac);

  ASSERT_NEAR(update.lambda[0], 0.5, eps<double>);
  ASSERT_NEAR(update.value, sqrt2, eps<double>);
}

TEST (updates_tri, mp0_basic_factoring_test) {
  double u0, u1, s, s0, s1, h, s_fac;
  vec<double, 2> p0, p1, p_fac;

  u0 = u1 = s = s0 = s1 = h = s_fac = 1;

  p0[0] = 1;
  p0[1] = 0;

  p1[0] = 0;
  p1[1] = 1;

  p_fac[0] = 1;
  p_fac[1] = 1;

  auto update = tri<MP0, 2>()(p0, p1, u0, u1, s, s0, s1, h, p_fac, s_fac);

  ASSERT_NEAR(update.lambda[0], 0.5, eps<double>);
  ASSERT_NEAR(update.value, sqrt2, eps<double>);
}

TEST (updates_tri, mp1_basic_factoring_test) {
  double u0, u1, s, s0, s1, h, s_fac;
  vec<double, 2> p0, p1, p_fac;

  u0 = u1 = s = s0 = s1 = h = s_fac = 1;

  p0[0] = 1;
  p0[1] = 0;

  p1[0] = 0;
  p1[1] = 1;

  p_fac[0] = 1;
  p_fac[1] = 1;

  auto update = tri<MP1, 2>()(p0, p1, u0, u1, s, s0, s1, h, p_fac, s_fac);

  ASSERT_NEAR(update.lambda[0], 0.5, eps<double>);
  ASSERT_NEAR(update.value, sqrt2, eps<double>);
}

TEST (updates_tri, factored_mp0_with_constant_slowness_works) {
  double u0, u1, s, s0, s1, h, s_fac;
  vec<double, 2> p0, p1, p_fac;

  u0 = u1 = s = s0 = s1 = h = s_fac = 1;

  p0[0] = p_fac[0] = 0;
  p0[1] = p_fac[1] = 1;

  p1[0] = 1;
  p1[1] = 1;

  double ufac = tri<MP0, 2>()(p0, p1, u0, u1, s, s0, s1, h, p_fac, s_fac).value;
  ASSERT_DOUBLE_EQ(ufac, tri12<MP0>(u0, u1, s, s0, s1, h).value);
}

TEST (updates_tri, factored_mp1_with_constant_slowness_works) {
  double u0, u1, s, s0, s1, h, s_fac;
  vec<double, 2> p0, p1, p_fac;

  u0 = u1 = s = s0 = s1 = h = s_fac = 1;

  p0[0] = p_fac[0] = 0;
  p0[1] = p_fac[1] = 1;

  p1[0] = 1;
  p1[1] = 1;

  double u = tri12<MP1>(u0, u1, s, s0, s1, h).value;
  double ufac = tri<MP1, 2>()(p0, p1, u0, u1, s, s0, s1, h, p_fac, s_fac).value;
  ASSERT_DOUBLE_EQ(u, ufac);
}

TEST (updates_tri, factored_rhr_with_constant_slowness_works) {
  double u, u_fac, u0, u1, s, s0, s1, h, s_fac;
  vec<double, 2> p0, p1, p_fac;

  {
    u0 = u1 = s = s0 = s1 = h = s_fac = 1;

    p0[0] = p_fac[0] = 0;
    p0[1] = p_fac[1] = 1;

    p1[0] = 1;
    p1[1] = 1;

    u = tri12<RHR>(u0, u1, s, s0, s1, h).value;
    u_fac = tri<RHR, 2>()(p0, p1, u0, u1, s, s0, s1, h, p_fac, s_fac).value;
    ASSERT_DOUBLE_EQ(u, u_fac);
  }
}

TEST (updates_tri, rhr_non_bv_tri_update_works) {
  {
    vec<double, 2> p0 = {0, -1};
    vec<double, 2> p1 = {-1, -1};
    double U0 = 1, U1 = 0, s = 1, s0 = 1, s1 = 1, h = 1;
    double u0 = tri12<RHR>(U0, U1, s, s0, s1, h).value;
    double u1 = tri<RHR, 2>()(p0, p1, U0, U1, s, s0, s1, h).value;
    ASSERT_DOUBLE_EQ(u0, u1);
    ASSERT_DOUBLE_EQ(u0, sqrt2);
    ASSERT_DOUBLE_EQ(u1, sqrt2);
  }
}
