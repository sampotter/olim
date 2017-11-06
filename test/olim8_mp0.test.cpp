#include <cmath>

#include "olim8.hpp"
#include "speed_funcs.hpp"
#include "test.hpp"

void trivial_case_works() {
  olim8_mp0 m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 0.0);
}

void adjacent_update_works() {
  olim8_mp0 m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(1, 0), 0.5);
}

void neighboring_values_are_correct() {
  olim8_mp0 m {3, 3, 1};
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

void origin_test() {
  int M = 5, N = 5;
  double h = 0.1, x0 = h*(N - 1)/2., y0 = h*(M - 1)/2.;
  olim8_mp0 m {M, N, h, default_speed_func, x0, y0};
  m.add_boundary_node(2, 2);
  m.run();
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      IS_APPROX_EQUAL(
        m.get_value(i, j), default_speed_func_soln(h*j - x0, h*i - y0), 4e-2);
    }
  }
}

void sf1_single_row_test() {
  int N = 1001;
  double h = 1.0/(N - 1);
  olim8_mp0 m {1, N, h, (speed_func) s1};
  m.add_boundary_node(0, 0);
  m.run();
  for (int j = N - 10; j < N; ++j) {
    double U = m.get_value(0, j);
    double u = f1(h*j, 0);
    printf("%g\n", fabs(u - U)/fabs(u));
    IS_APPROX_EQUAL(u, U, 1e-2);
  }
}

void sf1_test() {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp0 m {M, N, h, (speed_func) s1};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f1(h*j, h*i);
      printf("%g\n", fabs(u - U)/fabs(u));
      IS_APPROX_EQUAL(u, U, 1e-1);
    }
  }
}

int main() {
  trivial_case_works();
  adjacent_update_works();
  neighboring_values_are_correct();
  origin_test();
  sf1_single_row_test();
  sf1_test();
}
