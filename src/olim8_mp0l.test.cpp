#include <cmath>

#include "olim8_mp0c.hpp"
#include "olim8_mp0l.hpp"
#include "olim8_rhr.hpp"
#include "speed_funcs.hpp"
#include "test.hpp"

void trivial_case_works() {
  olim8_mp0l m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  test::is_approx_equal(m.get_value(0, 0), 0.0);
}

void adjacent_update_works() {
  olim8_mp0l m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  test::is_approx_equal(m.get_value(1, 0), 0.5);
}

void neighboring_values_are_correct() {
  olim8_mp0l m {3, 3, 1};
  m.add_boundary_node(1, 1);
  m.run();
  double gt[] = {
    std::sqrt(2),
    1,
    std::sqrt(2),
    1,
    0,
    1,
    std::sqrt(2),
    1,
    std::sqrt(2)
  };
  int k = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      test::is_approx_equal(gt[k++], m.get_value(i, j), 1e-15);
    }
  }
}

void degenerates_to_olim8_rhr_on_constant_speed_func() {
  int n = 5;
  double h = 2.0/(n - 1);

  olim8_mp0l m1 {n, n, h, default_speed_func, 1, 1};
  m1.add_boundary_node((n - 1)/2, (n - 1)/2);
  m1.run();

  olim8_rhr m2 {n, n, h, default_speed_func, 1, 1};
  m2.add_boundary_node((n - 1)/2, (n - 1)/2);
  m2.run();

  double x, y;
  for (int i = 0; i < n; ++i) {
    y = h*i - 1;
    for (int j = 0; j < n; ++j) {
      x = h*j - 1;
      test::is_approx_equal(m1.get_value(i, j), m2.get_value(i, j), 1e-15);
    }
  }
}

void piecewise_linear_test() {
  int n = 5;
  double h = 2.0/(n - 1);
  
  olim8_mp0l m1 {n, n, h, s2, 1, 1};
  m1.add_boundary_node((n - 1)/2, (n - 1)/2);
  m1.run();
  
  olim8_mp0c m2 {n, n, h, s2, 1, 1};
  m2.add_boundary_node((n - 1)/2, (n - 1)/2);
  m2.run();
  
  double x, y;
  for (int i = 0; i < n; ++i) {
    y = h*i - 1;
    for (int j = 0; j < n; ++j) {
      x = h*j - 1;
      test::is_approx_equal(m1.get_value(i, j), m2.get_value(i, j));
      test::is_approx_equal(m1.get_value(i, j), f2(x, y), 1e-15);
    }
  }
}

int main() {
  // trivial_case_works();
  // adjacent_update_works();
  // neighboring_values_are_correct();
  // degenerates_to_olim8_rhr_on_constant_speed_func();
  piecewise_linear_test();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
