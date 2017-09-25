#include "common.defs.hpp"
#include "conics.hpp"
#include "test.hpp"
#include "olim_update_rules.hpp"

#include <cmath>

static olim3d_rhr_update_rules<arma_rootfinder> updates;

void line1_works() {
  IS_APPROX_EQUAL(updates.line1(0.0, 1.0, 1.0, 1.0), 1.0);
}

void line2_works() {
  IS_APPROX_EQUAL(updates.line2(0.0, 1.0, 1.0, 1.0), sqrt(2));
}

void line3_works() {
  IS_APPROX_EQUAL(updates.line3(0.0, 1.0, 1.0, 1.0), sqrt(3));
}

void tri11_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, h = 1;
  double val01 = updates.tri11(u0, u1, s, s, s, h);
  double val10 = updates.tri11(u1, u0, s, s, s, h);
  IS_APPROX_EQUAL(val01, val10);
}

void tri11_works() {
  IS_APPROX_EQUAL(updates.tri11(0.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt(2)/2);
  IS_APPROX_EQUAL(updates.tri11(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
  IS_APPROX_EQUAL(updates.tri11(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

void tri12_works() {
  IS_APPROX_EQUAL(updates.tri12(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt(2));
  IS_APPROX_EQUAL(updates.tri12(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

void tri13_works() {
  double val = 1 - cos(pi/4) + sqrt(1 + 2*cos(pi/4)*cos(pi/4));
  IS_APPROX_EQUAL(updates.tri13(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), val);
  IS_APPROX_EQUAL(updates.tri13(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0);
}

void tri22_works() {
  IS_APPROX_EQUAL(updates.tri22(0.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt(6)/2);
  IS_APPROX_EQUAL(updates.tri22(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), sqrt(2));
  IS_APPROX_EQUAL(updates.tri22(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt(2));
}

void tri22_is_symmetric() {
  double u0 = 0.0, u1 = 0.1, s = 1, h = 1;
  double val01 = updates.tri22(u0, u1, s, s, s, h);
  double val10 = updates.tri22(u1, u0, s, s, s, h);
  IS_APPROX_EQUAL(val01, val10);
}

void tri23_works() {
  IS_APPROX_EQUAL(updates.tri23(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt(3));
  IS_APPROX_EQUAL(updates.tri23(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), sqrt(2));
}

void tetra111_is_symmetric() {
  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1, s1 = 1, s2 = 1, h = 1;
  double val012 = updates.tetra111(u0, u1, u2, s, s0, s1, s2, h);
  double val120 = updates.tetra111(u1, u2, u0, s, s1, s2, s0, h);
  double val201 = updates.tetra111(u2, u0, u1, s, s2, s0, s1, h);
  double val210 = updates.tetra111(u2, u1, u0, s, s2, s1, s0, h);
  double val102 = updates.tetra111(u1, u0, u2, s, s1, s0, s2, h);
  double val021 = updates.tetra111(u0, u2, u1, s, s0, s2, s1, h);
  IS_APPROX_EQUAL(val012, val120);
  IS_APPROX_EQUAL(val120, val201);
  IS_APPROX_EQUAL(val201, val210);
  IS_APPROX_EQUAL(val210, val102);
  IS_APPROX_EQUAL(val102, val021);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = updates.tetra111(u0, u1, u2, s, s0, s1, s2, h);
  val120 = updates.tetra111(u1, u2, u0, s, s1, s2, s0, h);
  val201 = updates.tetra111(u2, u0, u1, s, s2, s0, s1, h);
  val210 = updates.tetra111(u2, u1, u0, s, s2, s1, s0, h);
  val102 = updates.tetra111(u1, u0, u2, s, s1, s0, s2, h);
  val021 = updates.tetra111(u0, u2, u1, s, s0, s2, s1, h);
  IS_APPROX_EQUAL(val012, val120);
  IS_APPROX_EQUAL(val120, val201);
  IS_APPROX_EQUAL(val201, val210);
  IS_APPROX_EQUAL(val210, val102);
  IS_APPROX_EQUAL(val102, val021);
}

void tetra111_works() {
  IS_APPROX_EQUAL(updates.tetra111(1, 0, 0, 1, 1, 1, 1, 1), 1.0);
  IS_APPROX_EQUAL(updates.tetra111(0, 1, 0, 1, 1, 1, 1, 1), 1.0);
  IS_APPROX_EQUAL(updates.tetra111(0, 0, 1, 1, 1, 1, 1, 1), 1.0);
  IS_APPROX_EQUAL(updates.tetra111(0, 0, 0, 1, 1, 1, 1, 1), sqrt(3)/3);
}

void tetra122_is_symmetric() {
  double u0 = 2.4;
  double u1 = 2.2;
  double u2 = 2.0;
  double s = 1, s0 = 1, s1 = 1, s2 = 1, h = 1;
  double val012 = updates.tetra122(u0, u1, u2, s, s0, s1, s2, h);
  double val021 = updates.tetra122(u0, u2, u1, s, s0, s2, s1, h);
  IS_APPROX_EQUAL(val012, val021);
}

void tetra222_is_symmetric() {
  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1, s1 = 1, s2 = 1, h = 1;
  double val012 = updates.tetra222(u0, u1, u2, s, s0, s1, s2, h);
  double val120 = updates.tetra222(u1, u2, u0, s, s1, s2, s0, h);
  double val201 = updates.tetra222(u2, u0, u1, s, s2, s0, s1, h);
  double val210 = updates.tetra222(u2, u1, u0, s, s2, s1, s0, h);
  double val102 = updates.tetra222(u1, u0, u2, s, s1, s0, s2, h);
  double val021 = updates.tetra222(u0, u2, u1, s, s0, s2, s1, h);
  IS_APPROX_EQUAL(val012, val120);
  IS_APPROX_EQUAL(val120, val201);
  IS_APPROX_EQUAL(val201, val210);
  IS_APPROX_EQUAL(val210, val102);
  IS_APPROX_EQUAL(val102, val021);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = updates.tetra222(u0, u1, u2, s, s0, s1, s2, h);
  val120 = updates.tetra222(u1, u2, u0, s, s1, s2, s0, h);
  val201 = updates.tetra222(u2, u0, u1, s, s2, s0, s1, h);
  val210 = updates.tetra222(u2, u1, u0, s, s2, s1, s0, h);
  val102 = updates.tetra222(u1, u0, u2, s, s1, s0, s2, h);
  val021 = updates.tetra222(u0, u2, u1, s, s0, s2, s1, h);
  IS_APPROX_EQUAL(val012, val120);
  IS_APPROX_EQUAL(val120, val201);
  IS_APPROX_EQUAL(val201, val210);
  IS_APPROX_EQUAL(val210, val102);
  IS_APPROX_EQUAL(val102, val021);
}

void tetra222_works() {
  IS_APPROX_EQUAL(updates.tetra222(0, 0, 0, 1, 1, 1, 1, 1), 2.0/sqrt(3));
}

int main() {
  line1_works();
  line2_works();
  line3_works();
  tri11_is_symmetric();
  tri11_works();
  tri12_works();
  tri13_works();
  tri22_is_symmetric();
  tri22_works();
  tri23_works();
  tetra111_is_symmetric();
  tetra111_works();
  tetra122_is_symmetric();
  tetra222_is_symmetric();
  tetra222_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
