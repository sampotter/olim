#include "common.defs.hpp"
#include "common.macros.hpp"
#include "test.hpp"
#include "update_rules.tri_updates.hpp"

static update_rules::rhr_tri_updates<false> updates;

void tri11_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, s0 = 1, s1 = 1, h = 1;
  double val01 = updates.tri11(u0, u1, s, s0, s1, h);
  double val10 = updates.tri11(u1, u0, s, s1, s0, h);
  IS_APPROX_EQUAL(val01, val10);
}

void tri11_works() {
  IS_APPROX_EQUAL(updates.tri11(0.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt2/2);
  IS_APPROX_EQUAL(updates.tri11(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
  IS_APPROX_EQUAL(updates.tri11(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

void tri12_works() {
  // TODO: add interior point test
  IS_APPROX_EQUAL(updates.tri12(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt2);
  IS_APPROX_EQUAL(updates.tri12(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

void tri13_works() {
  // TODO: add interior point test
  IS_APPROX_EQUAL(updates.tri13(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), INF(double));
  IS_APPROX_EQUAL(updates.tri13(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), INF(double));
}

void tri22_works() {
  IS_APPROX_EQUAL(updates.tri22(0.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt6/2);
  IS_APPROX_EQUAL(updates.tri22(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), INF(double));
  IS_APPROX_EQUAL(updates.tri22(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), INF(double));
}

void tri22_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, s0 = 1.0, s1 = 1.0, h = 1;
  double val01 = updates.tri22(u0, u1, s, s0, s1, h);
  double val10 = updates.tri22(u1, u0, s, s1, s0, h);
  IS_APPROX_EQUAL(val01, val10);
}

void tri23_works() {
  IS_APPROX_EQUAL(updates.tri23(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), INF(double));
  IS_APPROX_EQUAL(updates.tri23(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), INF(double));
}

int main() {
  tri11_is_symmetric();
  tri11_works();
  tri12_works();
  tri13_works();
  tri22_is_symmetric();
  tri22_works();
  tri23_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
