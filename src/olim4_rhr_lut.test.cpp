#include "olim4_rhr_lut.hpp"
#include "speed_funcs.hpp"
#include "test.hpp"

void trivial_case_works() {
  olim4_rhr_lut m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  test::is_approx_equal(m.get_value(0, 0), 0.0);
}

void adjacent_update_works() {
  olim4_rhr_lut m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  test::is_approx_equal(m.get_value(1, 0), 0.5);
}

void neighboring_values_are_correct() {
  olim4_rhr_lut m {3, 3, 1};
  m.add_boundary_node(1, 1);
  m.run();
  double gt[] = {
    1.707106781186547,
    1,
    1.707106781186547,
    1,
    0,
    1,
    1.707106781186547,
    1,
    1.707106781186547
  };
  int k = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      test::is_approx_equal(gt[k++], m.get_value(i, j));
    }
  }
}

void masha_s_values_are_correct() {
  double gt[] = {
    0.02939463262323,
    0.01540058973304,
    0.02037282789484,
    0.01089789836732,
    2.030579675453e-07,
    0.01070992017687,
    0.02286521936613,
    0.01465429342172,
    0.02703978500521,
  };

  olim4_rhr_lut m {21, 21, 0.1, s7, 1, 1};
  m.add_boundary_node(10, 10, 2.030579675453e-07);
  m.run();

  int k = 0;
  for (int i = 9; i <= 11; ++i) {
    for (int j = 9; j <= 11; ++j) {
      test::is_approx_equal(gt[k++], m.get_value(i, j));
    }
  }
}

void masha_small_test() {
  double gt[] = {
    0.02939463262323,
    0.01540058973304,
    0.02037282789484,
    0.01089789836732,
    2.030579675453e-07,
    0.01070992017687,
    0.02286521936613,
    0.01465429342172,
    0.02703978500521,
  };
  double h = 0.1, x0 = 0.1, y0 = 0.1;
  olim4_rhr_lut m {3, 3, h, s7, x0, y0};
  m.add_boundary_node(1, 1, 2.030579675453e-07);
  m.run();
  for (int i = 0, k = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      test::is_approx_equal(m.get_value(i, j), gt[k++]);
    }
  }
}

void rectangular_domain_works() {
  double h = 0.1, x0 = 1, y0 = 1;
  olim4_rhr_lut m {11, 21, h, default_speed_func, x0, y0};
  m.add_boundary_node(5, 10);
  m.run();
  test::is_approx_equal(m.get_value(0, 0), 1.17825);
}

int main() {
  trivial_case_works();
  adjacent_update_works();
  neighboring_values_are_correct();
  masha_s_values_are_correct();
  masha_small_test();
  rectangular_domain_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
