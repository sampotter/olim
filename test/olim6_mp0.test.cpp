#include "test.hpp"

#include "basic_marcher_3d.hpp"
#include "olim4.hpp"
#include "olim6.hpp"

void planes_are_correct() {
  int n = 11;
  double h = 1.0/(n/2);
  
  olim4_mp0 m4 {n, n, h, default_speed_func, 1, 1};
  m4.add_boundary_node(n/2, n/2);
  m4.run();
  
  olim6_mp0 m6 {n, n, n, h, default_speed_func_3d, 1, 1, 1};
  m6.add_boundary_node(n/2, n/2, n/2);
  m6.run();

  // Check that planes are correct:
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double m4value = m4.get_value(i, j);
      IS_APPROX_EQUAL(m4value, m6.get_value(i, j, n/2));
      IS_APPROX_EQUAL(m4value, m6.get_value(i, n/2, j));
      IS_APPROX_EQUAL(m4value, m6.get_value(n/2, i, j));
    }
  }
}

void result_is_symmetric() {
  int n = 5;
  olim6_mp0 m {n, n, n, 1.0, default_speed_func_3d, 1.0, 1.0, 1.0};
  m.add_boundary_node(n/2, n/2, n/2);
  m.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0, k_ = n - 1; k < n; ++k, --k_) {
        IS_APPROX_EQUAL(m.get_value(i, j, k), m.get_value(i, j, k_));
      }
    }
  }

  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      for (int j = 0, j_ = n - 1; j < n; ++j, --j_) {
        IS_APPROX_EQUAL(m.get_value(i, j, k), m.get_value(i, j_, k));
      }
    }
  }

  for (int j = 0; j < n; ++j) {
    for (int k = 0; k < n; ++k) {
      for (int i = 0, i_ = n - 1; i < n; ++i, --i_) {
        IS_APPROX_EQUAL(m.get_value(i, j, k), m.get_value(i_, j, k));
      }
    }
  }
}

int main() {
  planes_are_correct();
  result_is_symmetric();
}
