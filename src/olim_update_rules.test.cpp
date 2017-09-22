#include "conics.hpp"
#include "test.hpp"
#include "olim_update_rules.hpp"

#include <cmath>

constexpr double pi = 3.141592653589793;

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

void tri23_works() {
  IS_APPROX_EQUAL(updates.tri23(1.0, 0.0, 1.0, 1.0, 1.0, 1.0), sqrt(3));
  IS_APPROX_EQUAL(updates.tri23(0.0, 1.0, 1.0, 1.0, 1.0, 1.0), sqrt(2));
}

int main() {
  line1_works();
  line2_works();
  line3_works();
  tri11_works();
  tri12_works();
  tri13_works();
  tri23_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
