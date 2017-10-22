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

void slightly_more_involved() {
  olim8_rhr m {2, 2, 1};
  m.add_boundary_node(0, 0);
  m.add_boundary_node(1, 0);
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

void run_on_full_neighborhood() {
  olim8_rhr m {3, 3, 1};
  m.add_boundary_node(0, 0);
  m.add_boundary_node(0, 1);
  m.add_boundary_node(0, 2);
  m.add_boundary_node(1, 0);
  m.add_boundary_node(1, 2);
  m.add_boundary_node(2, 0);
  m.add_boundary_node(2, 1);
  m.add_boundary_node(2, 2);
  m.run();
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

void s1_single_row_test() {
  int N = 1001;
  double h = 1.0/(N - 1);
  olim8_rhr m {1, N, h, (speed_func) s1};
  m.add_boundary_node(0, 0);
  m.run();
  double maxerr = 0;
  for (int j = 0; j < N; ++j) {
    double U = m.get_value(0, j);
    double u = f1(h*j, 0);
    maxerr = std::max(maxerr, std::fabs(U - u)/std::fabs(u));
  }
  std::cout << maxerr << std::endl;
}

void s1_test1() {
  int M = 1001, N = M;
  double h = 1.0/(M - 1);
  olim8_rhr m {M, N, h, (speed_func) s1};
  m.add_boundary_node(0, 0);
  m.run();
  double maxerr = 0;
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f1(h*j, h*i);
      if (u != 0) {
        maxerr = std::max(maxerr, std::fabs(U - u)/std::fabs(u));
      }
    }
  }
  std::cout << maxerr << std::endl;
}

void s1_test2() {
  int M = 201, N = M;
  double h = 2.0/(M - 1);
  olim8_rhr m {M, N, h, (speed_func) s1};
  m.add_boundary_node(100, 100);
  m.run();
  double maxerr = 0;
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      double U = m.get_value(i, j);
      double u = f1(h*j, h*i);
      if (u != 0) {
        maxerr = std::max(maxerr, std::fabs(U - u)/std::fabs(u));
      }
    }
  }
  std::cout << maxerr << std::endl;
}

int main() {
  trivial_case_works();
  adjacent_update_works();
  neighboring_values_are_correct();
  slightly_more_involved();
  slightly_more_involved_2();
  run_on_full_neighborhood();
  maria_test();
  origin_test();
  s1_single_row_test();
  s1_test1();
  s1_test2();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
