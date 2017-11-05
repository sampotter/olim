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
  {
    u0 = 0.1, u1 = 0, s = 1, s0 = 1.2, s1 = 1.1, h = 0.9, T = 0.73102427;
    IS_APPROX_EQUAL(mp1c.tri11(u0, u1, s, s0, s1, h), T);
    IS_APPROX_EQUAL(mp1u.tri11(u0, u1, s, s0, s1, h), T);
  }
  {
    u0 = 0, u1 = 0.1, s = 1, s0 = 1.1, s1 = 1.3, h = 1.2, T = 0.978772;
    IS_APPROX_EQUAL(mp1c.tri11(u0, u1, s, s0, s1, h), T, 1e-3);
    IS_APPROX_EQUAL(mp1u.tri11(u0, u1, s, s0, s1, h), T, 1e-3);
  }
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
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
