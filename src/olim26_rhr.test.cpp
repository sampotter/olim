#include "test.hpp"

#include <cmath>

#include "olim8_rhr.hpp"
#include "olim26_rhr.hpp"

constexpr double pi = 3.141592653589793;

static olim26_rhr_update_rules<arma_rootfinder> updates;

void line1_works() {
  IS_APPROX_EQUAL(updates.line1(0.0, 1.0, 1.0, 1.0), 1.0);
}

void line2_works() {
  IS_APPROX_EQUAL(updates.line2(0.0, 1.0, 1.0, 1.0), sqrt(2));
}

void line3_works() {
  IS_APPROX_EQUAL(updates.line3(0.0, 1.0, 1.0, 1.0), sqrt(3));
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

void planes_are_correct() {
  int n = 21;
  double h = 1.0/(n/2);
  
  olim8_rhr m8 {n, n, h, default_speed_func, 1, 1};
  m8.add_boundary_node(n/2, n/2);
  m8.run();
  
  olim26_rhr_arma m26 {n, n, n, h, default_speed_func_3d, 1, 1, 1};
  m26.add_boundary_node(n/2, n/2, n/2);
  m26.run();

  // Check that planes are correct:
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double m8value = m8.get_value(i, j);
      IS_APPROX_EQUAL(m8value, m26.get_value(i, j, n/2));
      IS_APPROX_EQUAL(m8value, m26.get_value(i, n/2, j));
      IS_APPROX_EQUAL(m8value, m26.get_value(n/2, i, j));
    }
  }
}

int main() {
  line1_works();
  line2_works();
  line3_works();
  tri12_works();
  tri13_works();
  tri23_works();
  planes_are_correct();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
