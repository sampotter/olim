#include "common.defs.hpp"
#include "common.macros.hpp"
#include "test.hpp"
#include "update_rules.tri_updates.hpp"

update_rules::rhr_tri_updates<true> rhrc; // (c)onstrained
update_rules::rhr_tri_updates<false> rhru; // (u)nconstrained

void rhr_tri11_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, s0 = 1, s1 = 1, h = 1;
  IS_APPROX_EQUAL(
    rhrc.tri11(u0, u1, s, s0, s1, h),
    rhrc.tri11(u1, u0, s, s1, s0, h));
  IS_APPROX_EQUAL(
    rhru.tri11(u0, u1, s, s0, s1, h),
    rhru.tri11(u1, u0, s, s1, s0, h));
}

void rhr_tri11_works() {
  IS_APPROX_EQUAL(rhrc.tri11(0.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt2/2);
  IS_APPROX_EQUAL(rhrc.tri11(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
  IS_APPROX_EQUAL(rhrc.tri11(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

void rhr_tri12_works() {
  // TODO: add interior point test
  IS_APPROX_EQUAL(rhrc.tri12(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt2);
  IS_APPROX_EQUAL(rhrc.tri12(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);

  double u0, u1, s, s0, s1, h;

  u0 = 0.923028, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri12(u0, u1, s, s0, s1, h), INF(double));

  u0 = 0.852848, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri12(u0, u1, s, s0, s1, h), 1.05597249);

  u0 = 0.923028, u1 = 0.817579, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri12(u0, u1, s, s0, s1, h), 1.14970069);

  u0 = 0.883639, u1 = 0.817579, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri12(u0, u1, s, s0, s1, h), 1.12475323);

  u0 = 0.852848, u1 = 0.65974, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri12(u0, u1, s, s0, s1, h), INF(double));

  u0 = 0.883639, u1 = 0.65974, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri12(u0, u1, s, s0, s1, h), INF(double));
}

void rhr_tri13_works() {
  // TODO: add interior point test
  IS_APPROX_EQUAL(rhru.tri13(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), INF(double));
  IS_APPROX_EQUAL(rhru.tri13(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), INF(double));
}

void rhr_tri22_works() {
  IS_APPROX_EQUAL(rhru.tri22(0.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt6/2);
  IS_APPROX_EQUAL(rhru.tri22(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), INF(double));
  IS_APPROX_EQUAL(rhru.tri22(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), INF(double));

  double u0, u1, s, s0, s1, h;

  u0 = 0.65974, u1 = 0.817579, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri22(u0, u1, s, s0, s1, h), 1.012639677310786);

  u0 = 0.65974, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri22(u0, u1, s, s0, s1, h), 0.986849397845339);

  u0 = 0.707107, u1 = 0.817579, s = 1, s0 = 1, s1 = 1, h = 0.25;
  IS_APPROX_EQUAL(rhru.tri22(u0, u1, s, s0, s1, h), 1.053198553345643);
}

void rhr_tri22_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, s0 = 1.0, s1 = 1.0, h = 1;
  double val01 = rhru.tri22(u0, u1, s, s0, s1, h);
  double val10 = rhru.tri22(u1, u0, s, s1, s0, h);
  IS_APPROX_EQUAL(val01, val10);
}

void rhr_tri23_works() {
  IS_APPROX_EQUAL(rhru.tri23(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), INF(double));
  IS_APPROX_EQUAL(rhru.tri23(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), INF(double));
}

// TODO: test constrained rhr

update_rules::mp1_tri_updates<true> mp1c; // (c)onstrained
update_rules::mp1_tri_updates<false> mp1u; // (u)nconstrained

void mp1_tri11_works() {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.1, u1 = 0, s = 1, s0 = 1.2, s1 = 1.1, h = 0.9, T = 0.73093624;
  IS_APPROX_EQUAL(mp1c.tri11(u0, u1, s, s0, s1, h), T);
  IS_APPROX_EQUAL(mp1u.tri11(u0, u1, s, s0, s1, h), T);

  u0 = 0, u1 = 0.1, s = 1, s0 = 1.1, s1 = 1.3, h = 1.2, T = 0.97877243;
  IS_APPROX_EQUAL(mp1c.tri11(u0, u1, s, s0, s1, h), T);
  IS_APPROX_EQUAL(mp1u.tri11(u0, u1, s, s0, s1, h), T);

  u0 = 1.0, u1 = 0.853553, s = 1.0, s0 = 1.0, s1 = 1.0, h = 0.5, T = 1.27266;
  IS_APPROX_EQUAL(mp1c.tri11(u0, u1, s, s0, s1, h), T, 4e-6);
  IS_APPROX_EQUAL(mp1u.tri11(u0, u1, s, s0, s1, h), T, 4e-6);
}

void mp1_tri12_works() {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.46487865669951, u1 = 0.4, s = 1, s0 = 1, s1 = 1, h = 0.2;
  T = 0.6540631166978595;
  IS_APPROX_EQUAL(mp1u.tri12(u0, u1, s, s0, s1, h), T);
}

void mp1_tri13_works() {
  double u0, u1, s, s0, s1, h;

  u0 = 0, u1 = 1.41421, s = 1, s0 = 1, s1 = 1, h = 1;
  IS_APPROX_EQUAL(mp1u.tri13(u0, u1, s, s0, s1, h), INF(double));
}

void mp1_tri22_works() {
  double u0, u1, s, s0, s1, h, T;

  u0 = 0.65974, u1 = 0.707107, s = 1, s0 = 1, s1 = 1, h = 0.25;
  T = 0.986849397845339;
  IS_APPROX_EQUAL(mp1u.tri22(u0, u1, s, s0, s1, h), T);
}

int main() {
  rhr_tri11_is_symmetric();
  rhr_tri11_works();
  rhr_tri12_works();
  rhr_tri13_works();
  rhr_tri22_is_symmetric();
  rhr_tri22_works();
  rhr_tri23_works();
  mp1_tri11_works();
  mp1_tri12_works();
  mp1_tri13_works();
  mp1_tri22_works();
}
