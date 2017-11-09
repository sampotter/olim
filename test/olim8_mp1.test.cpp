#include "olim8.hpp"

#include "common.defs.hpp"
#include "olim.test.common.hpp"
#include "speed_funcs.hpp"
#include "test.hpp"

void neighboring_values_are_correct() {
  olim8_mp1 m {3, 3, 1};
  m.add_boundary_node(1, 1);
  m.run();
  double gt[] = {sqrt2, 1, sqrt2, 1, 0, 1, sqrt2, 1, sqrt2};
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
  olim8_mp1 m {M, N, h, default_speed_func, x0, y0};
  m.add_boundary_node(2, 2);
  m.run();
  for (int i = 0; i < 5; ++i) {
    for (int j = 0; j < 5; ++j) {
      IS_APPROX_EQUAL(
        m.get_value(i, j), default_speed_func_soln(h*j - x0, h*i - y0), 4e-2);
    }
  }
}

void s1_single_row_test() {
  int N = 1001;
  double h = 1.0/(N - 1);
  olim8_mp1 m {1, N, h, (speed_func) s1};
  m.add_boundary_node(0, 0);
  m.run();
  for (int j = N - 10; j < N; ++j) {
    double U = m.get_value(0, j);
    double u = f1(h*j, 0);
    IS_APPROX_EQUAL(u, U, 1e-2);
  }
}

void s1_test() {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, (speed_func) s1};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f1(h*j, h*i);
      IS_APPROX_EQUAL(u, U, 1e-1);
    }
  }
}

void s2_test() {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, s2};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f2(h*j, h*i);
      IS_APPROX_EQUAL(u, U, 1e-1);
    }
  }
}

void s3_test() {
  int M = 101, N = M;
  double h = 1.0/(M - 1);
  olim8_mp1 m {M, N, h, s3};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = M - 3; i < M; ++i) {
    for (int j = N - 3; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f3(h*j, h*i);
      IS_APPROX_EQUAL(u, U, 1e-1);
    }
  }
}

void s_test(speed_func s, speed_func f, int n = 11) {
  double h = 1.0/(n - 1);
  olim8_mp1 m {n, n, h, s};
  m.add_boundary_node(0, 0);
  m.run();
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double U = m.get_value(i, j);
      double u = f(h*j, h*i);
      IS_APPROX_EQUAL(u, U, 1e-1);
    }
  }
}

int main() {
  using olim = olim8_mp1;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
  neighboring_values_are_correct();
  origin_test();
  s1_single_row_test();
  s1_test();
  s2_test();
  s3_test();
  s_test(s7, f7);
}
