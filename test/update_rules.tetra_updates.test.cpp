#include "common.defs.hpp"
#include "test.hpp"
#include "update_rules.tetra_updates.hpp"

#include <cmath>

static update_rules::rhr_tetra_updates updates;

void tetra111_is_symmetric() {
  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1, s1 = 1, s2 = 1;
  double h = 1;
  double val012 = updates.tetra111(u0, u1, u2, s, s0, s1, s2, h);
  double val120 = updates.tetra111(u1, u2, u0, s, s0, s1, s2, h);
  double val201 = updates.tetra111(u2, u0, u1, s, s0, s1, s2, h);
  double val210 = updates.tetra111(u2, u1, u0, s, s0, s1, s2, h);
  double val102 = updates.tetra111(u1, u0, u2, s, s0, s1, s2, h);
  double val021 = updates.tetra111(u0, u2, u1, s, s0, s1, s2, h);
  IS_APPROX_EQUAL(val012, val120);
  IS_APPROX_EQUAL(val120, val201);
  IS_APPROX_EQUAL(val201, val210);
  IS_APPROX_EQUAL(val210, val102);
  IS_APPROX_EQUAL(val102, val021);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = updates.tetra111(u0, u1, u2, s, s0, s1, s2, h);
  val120 = updates.tetra111(u1, u2, u0, s, s0, s1, s2, h);
  val201 = updates.tetra111(u2, u0, u1, s, s0, s1, s2, h);
  val210 = updates.tetra111(u2, u1, u0, s, s0, s1, s2, h);
  val102 = updates.tetra111(u1, u0, u2, s, s0, s1, s2, h);
  val021 = updates.tetra111(u0, u2, u1, s, s0, s1, s2, h);
  IS_APPROX_EQUAL(val012, val120);
  IS_APPROX_EQUAL(val120, val201);
  IS_APPROX_EQUAL(val201, val210);
  IS_APPROX_EQUAL(val210, val102);
  IS_APPROX_EQUAL(val102, val021);
}

void tetra111_works() {
  {
    double u0 = 1.867422661146497;
    double u1 = 1.872030476918322;
    double u2 = 1.874601455103048;
    double s = 1.0, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 1.0/7.0;
    double uhat = 1.95377665722661;
    IS_APPROX_EQUAL(updates.tetra111(u0, u1, u2, s, s0, s1, s2, h), uhat);
  }
  {
    double u0 = 0.6714002359494359;
    double u1 = 0.667974148883546;
    double u2 = 0.652837863557498;
    double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 0.02040816326530612;
    double uhat = 0.6726606175825081;
    IS_APPROX_EQUAL(updates.tetra111(u0, u1, u2, s, s0, s1, s2, h), uhat);
  }
  {
    double u0 = 0.8701508258299168;
    double u1 = 0.9034596879080383;
    double u2 = 0.8999244233472412;
    double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0;
    double h = 0.03703703703703703;
    double uhat = 0.9064782064785435;
    IS_APPROX_EQUAL(updates.tetra111(u0, u1, u2, s, s0, s1, s2, h), uhat);
  }
}

void tetra122_is_symmetric() {
  double u0 = 2.4;
  double u1 = 2.2;
  double u2 = 2.0;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = updates.tetra122(u0, u1, u2, s, s0, s1, s2, h);
  double val021 = updates.tetra122(u0, u2, u1, s, s0, s1, s2, h);
  IS_APPROX_EQUAL(val012, val021);
}

void tetra222_is_symmetric() {
  double u0 = 2.5453289254261224;
  double u1 = 2.5453289254261224;
  double u2 = 2.2844570503761732;
  double s = 1, s0 = 1.0, s1 = 1.0, s2 = 1.0, h = 1;
  double val012 = updates.tetra222(u0, u1, u2, s, s0, s1, s2, h);
  double val120 = updates.tetra222(u1, u2, u0, s, s0, s1, s2, h);
  double val201 = updates.tetra222(u2, u0, u1, s, s0, s1, s2, h);
  double val210 = updates.tetra222(u2, u1, u0, s, s0, s1, s2, h);
  double val102 = updates.tetra222(u1, u0, u2, s, s0, s1, s2, h);
  double val021 = updates.tetra222(u0, u2, u1, s, s0, s1, s2, h);
  IS_APPROX_EQUAL(val012, val120);
  IS_APPROX_EQUAL(val120, val201);
  IS_APPROX_EQUAL(val201, val210);
  IS_APPROX_EQUAL(val210, val102);
  IS_APPROX_EQUAL(val102, val021);

  u0 = 2.4;
  u1 = 2.2;
  u2 = 2.0;
  val012 = updates.tetra222(u0, u1, u2, s, s0, s1, s2, h);
  val120 = updates.tetra222(u1, u2, u0, s, s0, s1, s2, h);
  val201 = updates.tetra222(u2, u0, u1, s, s0, s1, s2, h);
  val210 = updates.tetra222(u2, u1, u0, s, s0, s1, s2, h);
  val102 = updates.tetra222(u1, u0, u2, s, s0, s1, s2, h);
  val021 = updates.tetra222(u0, u2, u1, s, s0, s1, s2, h);
  IS_APPROX_EQUAL(val012, val120);
  IS_APPROX_EQUAL(val120, val201);
  IS_APPROX_EQUAL(val201, val210);
  IS_APPROX_EQUAL(val210, val102);
  IS_APPROX_EQUAL(val102, val021);
}

void tetra222_works() {
  IS_APPROX_EQUAL(updates.tetra222(0, 0, 0, 1, 1, 1, 1, 1), 2.0/sqrt(3));
}

void tetra123_works() {
  IS_APPROX_EQUAL(updates.tetra123(0, 0, 0, 1, 1, 1, 1, 1), 1.0);
}

int main() {
  tetra111_is_symmetric();
  tetra111_works();
  tetra122_is_symmetric();
  tetra123_works();
  tetra222_is_symmetric();
  tetra222_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End: