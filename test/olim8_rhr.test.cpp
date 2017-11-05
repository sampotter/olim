#include "olim8.hpp"
#include "speed_funcs.hpp"
#include "test.hpp"

#include <cstdio>

void trivial_case_works() {
  olim8_rhr m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 0.0);
};

void adjacent_update_works() {
  olim8_rhr m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 0.0);
  IS_APPROX_EQUAL(m.get_value(1, 0), 0.5);
};

void neighboring_values_are_correct() {
  olim8_rhr m {3, 3, 1};
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
      IS_APPROX_EQUAL(gt[k++], m.get_value(i, j), 1e-15);
    }
  }
}

void maria_test() {
  olim8_rhr m {3, 2, 1};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 0.0);
  IS_APPROX_EQUAL(m.get_value(0, 1), 1.0);
  IS_APPROX_EQUAL(m.get_value(1, 0), 1.0);
  IS_APPROX_EQUAL(m.get_value(1, 1), std::sqrt(2.0));
  IS_APPROX_EQUAL(m.get_value(2, 0), 2.0);
  IS_APPROX_EQUAL(m.get_value(2, 1), 2.3243932834975496);
}

void origin_test() {
  int M = 5, N = 5;
  double h = 0.1, x0 = h*(N - 1)/2., y0 = h*(M - 1)/2.;
  olim8_rhr m {M, N, h, default_speed_func, x0, y0};
  m.add_boundary_node(2, 2);
  m.run();
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      IS_APPROX_EQUAL(
        m.get_value(i, j), default_speed_func_soln(h*j - x0, h*i - y0), 4e-2);
    }
  }
}

void slightly_more_involved() {
  olim8_rhr m {2, 2, 1};
  node nodes[2] = {{0, 0}, {1, 0}};
  m.add_boundary_nodes(nodes, 2);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 0.0);
  IS_APPROX_EQUAL(m.get_value(1, 0), 0.0);
  IS_APPROX_EQUAL(m.get_value(0, 1), 1.0);
  IS_APPROX_EQUAL(m.get_value(1, 1), 1.0);
}

void slightly_more_involved_2() {
  olim8_rhr m {2, 2, 1};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 0.0);
  IS_APPROX_EQUAL(m.get_value(1, 0), 1.0);
  IS_APPROX_EQUAL(m.get_value(0, 1), 1.0);
}

int main() {
  trivial_case_works();
  adjacent_update_works();
  neighboring_values_are_correct();
  maria_test();
  origin_test();
  slightly_more_involved();
  slightly_more_involved_2();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
