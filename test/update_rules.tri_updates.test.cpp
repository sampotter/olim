#include <gtest/gtest.h>

#include <src/config.hpp>

#include "common.defs.hpp"
#include "common.macros.hpp"
#include "update_rules.tri_updates.hpp"

using namespace update_rules;

mp0_tri_updates mp0;
mp1_tri_updates mp1;
rhr_tri_updates rhr;

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

#define TRI(tri_updates, p0, p1, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<p0> {}, ffvec<p1> {}).value

#define TRI11(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P01> {}, ffvec<P10> {}).value

#define TRI12(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P01> {}, ffvec<P11> {}).value

#define TRI13(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P001> {}, ffvec<P111> {}).value

#define TRI22(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P011> {}, ffvec<P101> {}).value

#define TRI23(tri_updates, u0, u1, s, s0, s1, h) \
  tri_updates.tri(u0, u1, s, s0, s1, h, ffvec<P011> {}, ffvec<P111> {}).value

/**
 * rhr tri11 tests
 */

TEST (tri_updates, rhr_tri11_basic_test) {
  double t = TRI11(rhr, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0);
  ASSERT_DOUBLE_EQ(t, sqrt2/2);
  t = TRI11(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0);
  ASSERT_DOUBLE_EQ(t, 1.0);
  t = TRI11(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0);
  ASSERT_DOUBLE_EQ(t, 1.0);
}

template <class rules>
void tri11_is_symmetric_with_constant_slowness() {
  rules r;
  double U0 = 0.0, U1 = 0.1, s = 1, s0 = 1, s1 = 1, h = 1, u, u1, u2;

  u = 0.7553367989832942;
  u1 = TRI11(r, U0, U1, s, s0, s1, h);
  u2 = TRI11(r, U1, U0, s, s1, s0, h);
  ASSERT_DOUBLE_EQ(u1, u);
  ASSERT_DOUBLE_EQ(u2, u);
  ASSERT_DOUBLE_EQ(u1, u2);
}

TEST (tri_updates, tri11_is_symmetric_with_constant_slowness) {
  tri11_is_symmetric_with_constant_slowness<mp0_tri_updates>();
  tri11_is_symmetric_with_constant_slowness<mp1_tri_updates>();
  tri11_is_symmetric_with_constant_slowness<rhr_tri_updates>();
}

TEST (tri_updates, tri11_rhr_is_symmetric_with_nonconstant_slowness) {
  rhr_tri_updates r;
  double U0 = 0.0, U1 = 0.1, s = 1, s0 = 1.2, s1 = 1.1, h = 1, u, u1, u2;
  u = 0.7553367989832942;
  u1 = TRI11(r, U0, U1, s, s0, s1, h);
  u2 = TRI11(r, U1, U0, s, s1, s0, h);
  ASSERT_DOUBLE_EQ(u1, u);
  ASSERT_DOUBLE_EQ(u2, u);
  ASSERT_DOUBLE_EQ(u1, u2);
}

TEST (tri_updates, tri11_mp0_is_symmetric_with_nonconstant_slowness) {
  mp0_tri_updates r;
  double U0 = 0.0, U1 = 0.1, s = 1, s0 = 1.2, s1 = 1.1, h = 1, u, u1, u2;
  u = 0.809661416377184;
  u1 = TRI11(r, U0, U1, s, s0, s1, h);
  u2 = TRI11(r, U1, U0, s, s1, s0, h);
  ASSERT_DOUBLE_EQ(u1, u);
  ASSERT_DOUBLE_EQ(u2, u);
  ASSERT_DOUBLE_EQ(u1, u2);
}

TEST (tri_updates, tri11_mp1_is_symmetric_with_nonconstant_slowness) {
  mp1_tri_updates r;
  double U0 = 0.0, U1 = 0.1, s = 1, s0 = 1.2, s1 = 1.1, h = 1, u, u1, u2;
  u = 0.809452960468600;
  u1 = TRI11(r, U0, U1, s, s0, s1, h);
  u2 = TRI11(r, U1, U0, s, s1, s0, h);
  ASSERT_DOUBLE_EQ(u1, u);
  ASSERT_DOUBLE_EQ(u2, u);
  ASSERT_DOUBLE_EQ(u1, u2);
}

TEST (tri_updates, rhr_tri12_works) {
  double t, u0, u1, s = 1, s0 = 1, s1 = 1, h = 0.25;

  t = TRI12(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0);
  ASSERT_DOUBLE_EQ(t, sqrt2);

  ASSERT_DOUBLE_EQ(TRI12(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);

  u0 = 0.923028, u1 = 0.707107;
  ASSERT_NEAR(TRI12(rhr, u0, u1, s, s0, s1, h), 1.060660390593274, 1e-14);

  u0 = 0.852848, u1 = 0.707107;
  ASSERT_NEAR(TRI12(rhr, u0, u1, s, s0, s1, h), 1.05597249, 6.2e-9);

  u0 = 0.923028, u1 = 0.817579;
  ASSERT_NEAR(TRI12(rhr, u0, u1, s, s0, s1, h), 1.14970069, 1.6e-11);

  u0 = 0.883639, u1 = 0.817579;
  ASSERT_NEAR(TRI12(rhr, u0, u1, s, s0, s1, h), 1.12475323, 9.32e-9);

  u0 = 0.852848, u1 = 0.65974;
  ASSERT_NEAR(TRI12(rhr, u0, u1, s, s0, s1, h), 1.013293390593274, 1e-14);

  u0 = 0.883639, u1 = 0.65974;
  ASSERT_NEAR(TRI12(rhr, u0, u1, s, s0, s1, h), 1.013293390593274, 1e-14);
}

TEST (tri_updates, rhr_tri13_works) {
  ASSERT_DOUBLE_EQ(TRI13(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.707106781186547);
  ASSERT_DOUBLE_EQ(TRI13(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

TEST (tri_updates, rhr_tri22_works) {
  ASSERT_DOUBLE_EQ(TRI22(rhr, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.224744871391589);
  ASSERT_DOUBLE_EQ(TRI22(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.414213562373095);
  ASSERT_DOUBLE_EQ(TRI22(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.414213562373095);

  double u0, u1, s = 1, s0 = 1, s1 = 1, h = 0.25;

  u0 = 0.65974, u1 = 0.817579;
  ASSERT_DOUBLE_EQ(TRI22(rhr, u0, u1, s, s0, s1, h), 1.012639677310786);

  u0 = 0.65974, u1 = 0.707107;
  ASSERT_DOUBLE_EQ(TRI22(rhr, u0, u1, s, s0, s1, h), 0.986849397845339);

  u0 = 0.707107, u1 = 0.817579;
  ASSERT_DOUBLE_EQ(TRI22(rhr, u0, u1, s, s0, s1, h), 1.053198553345643);
}

template <class rules>
void tri22_is_symmetric() {
  rules r;
  double u0 = 0.0, u1 = 0.1, s = 1, s0 = 1.0, s1 = 1.0, h = 1;
  double val01 = TRI22(r, u0, u1, s, s0, s1, h);
  double val10 = TRI22(r, u1, u0, s, s1, s0, h);
  ASSERT_DOUBLE_EQ(val01, val10);
}

TEST (tri_updates, tri22_is_symmetric_with_constant_slowness) {
  tri22_is_symmetric<mp0_tri_updates>();
  tri22_is_symmetric<mp1_tri_updates>();
  tri22_is_symmetric<rhr_tri_updates>();
}

TEST (tri_updates, rhr_tri23_works) {
  ASSERT_DOUBLE_EQ(TRI23(rhr, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.732050807568877);
  ASSERT_DOUBLE_EQ(TRI23(rhr, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.414213562373095);
}

TEST (tri_updates, mp1_tri11_works) {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.1, u1 = 0, s = 1, s0 = 1.2, s1 = 1.1, h = 0.9, T = 0.7309362364283433;
  ASSERT_NEAR(TRI11(mp1, u0, u1, s, s0, s1, h), T, 1e-13);
  ASSERT_NEAR(TRI11(mp1, u1, u0, s, s1, s0, h), T, 1e-13);

  u0 = 0, u1 = 0.1, s = 1, s0 = 1.1, s1 = 1.3, h = 1.2, T = 0.978772426584215;
  ASSERT_NEAR(TRI11(mp1, u0, u1, s, s0, s1, h), T, 1e-13);
  ASSERT_NEAR(TRI11(mp1, u1, u0, s, s1, s0, h), T, 1e-13);

  u0 = 1.0, u1 = 0.853553, s = 1.0, s0 = 1.0, s1 = 1.0, h = 0.5,
    T = 1.27266422607271;
  ASSERT_NEAR(TRI11(mp1, u0, u1, s, s0, s1, h), T, 1e-13);
  ASSERT_NEAR(TRI11(mp1, u1, u0, s, s1, s0, h), T, 1e-13);
}

TEST (tri_updates, mp1_tri12_works) {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.46487865669951, u1 = 0.4, s = 1, s0 = 1, s1 = 1, h = 0.2;
  T = 0.6540631166978595;
  ASSERT_DOUBLE_EQ(TRI12(mp1, u0, u1, s, s0, s1, h), T);
}

TEST (tri_updates, mp1_tri13_works) {
  double u0, u1, s, s0, s1, h;

  u0 = 2, u1 = 1, s = 1, s0 = 1, s1 = 1, h = 1;
  ASSERT_DOUBLE_EQ(TRI13(mp1, u0, u1, s, s0, s1, h), 2 + 1./std::sqrt(2));
}

TEST (tri_updates, mp1_tri22_works) {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.65974, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  T = 0.986849397845339;
  ASSERT_DOUBLE_EQ(TRI22(mp1, u0, u1, s, s0, s1, h), T);

  u0 = 0.59457611725233628;
  u1 = 0.56102682919361013;
  s = 0.0025050133959455545;
  s0 = 0.0025050133959455545;
  s1 = 0.2382400185837108;
  h = 0.5;
  T = 0.5963474292115725;
  ASSERT_DOUBLE_EQ(TRI22(mp1, u0, u1, s, s0, s1, h), T);
}

TEST (tri_updates, rhr_basic_factoring_test) {
  double u0, u1, s, s0, s1, h, p0[2], p1[2], p_fac[2], s_fac;

  u0 = u1 = s = s0 = s1 = h = s_fac = 1;

  p0[0] = 1;
  p0[1] = 0;

  p1[0] = 0;
  p1[1] = 1;

  p_fac[0] = 1;
  p_fac[1] = 1;

  auto update = rhr.tri<2>(u0, u1, s, s0, s1, h, p0, p1, p_fac, s_fac);

  ASSERT_NEAR(update.lambda[0], 0.5, EPS(double));
  ASSERT_NEAR(update.value, sqrt2, EPS(double));
}

TEST (tri_updates, factored_mp0_with_constant_slowness_works) {
  double u0, u1, s, s0, s1, h, p0[2], p1[2], p_fac[2], s_fac;

  u0 = u1 = s = s0 = s1 = h = s_fac = 1;

  p0[0] = p_fac[0] = 0;
  p0[1] = p_fac[1] = 1;

  p1[0] = 1;
  p1[1] = 1;

  ASSERT_DOUBLE_EQ(
    TRI12(mp0, u0, u1, s, s0, s1, h),
    mp0.tri<2>(u0, u1, s, s0, s1, h, p0, p1, p_fac, s_fac).value);
}

TEST (tri_updates, factored_mp1_with_constant_slowness_works) {
  double u0, u1, s, s0, s1, h, p0[2], p1[2], p_fac[2], s_fac;

  u0 = u1 = s = s0 = s1 = h = s_fac = 1;

  p0[0] = p_fac[0] = 0;
  p0[1] = p_fac[1] = 1;

  p1[0] = 1;
  p1[1] = 1;

  ASSERT_DOUBLE_EQ(
    TRI12(mp1, u0, u1, s, s0, s1, h),
    mp1.tri<2>(u0, u1, s, s0, s1, h, p0, p1, p_fac, s_fac).value);
}

TEST (tri_updates, factored_rhr_with_constant_slowness_works) {
  double u0, u1, s, s0, s1, h, p0[2], p1[2], p_fac[2], s_fac;

  {
    u0 = u1 = s = s0 = s1 = h = s_fac = 1;

    p0[0] = p_fac[0] = 0;
    p0[1] = p_fac[1] = 1;

    p1[0] = 1;
    p1[1] = 1;

    ASSERT_DOUBLE_EQ(
      TRI12(rhr, u0, u1, s, s0, s1, h),
      rhr.tri<2>(u0, u1, s, s0, s1, h, p0, p1, p_fac, s_fac).value);
  }

  {
    u0 = 0.44169958130222325;
    u1 = 0.44726430952545132;
    s = 0.34651473755893192;
    s0 = 0.39000936008666987;
    s1 = 0.37029340913370545;
    h = 0.040000000000000001;
    s_fac = 1;

    p0[0] = 1;
    p0[1] = 1;
    
    p1[0] = 0;
    p1[1] = 1;

    p_fac[0] = 11;
    p_fac[1] = 14;

    rhr.tri<2>(u0, u1, s, s0, s1, h, p0, p1, p_fac, s_fac);
  }
}
