#include "test.hpp"

#include "olim4_rhr.hpp"
#include "olim6_rhr.hpp"

void quadrants_are_correct() {
  int n = 2;
  double h = 1;

  /*
   * Tests for quadrants in the x-y plane:
   */
  {
    olim6_rhr_arma m {n, n, 1, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0 + sqrt(2)/2.0);
  }
  {
    olim6_rhr_arma m {n, n, 1, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim6_rhr_arma m {n, n, 1, h, default_speed_func_3d, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim6_rhr_arma m {n, n, 1, h, default_speed_func_3d, 1, 1, 0};
    m.add_boundary_node(1, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 0.0);
  }

  /**
   * Tests for quadrants in the x-z plane:
   */
  {
    olim6_rhr_arma m {n, 1, n, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0 + sqrt(2)/2.0);
  }
  {
    olim6_rhr_arma m {n, 1, n, h, default_speed_func_3d, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0);
  }
  {
    olim6_rhr_arma m {n, 1, n, h, default_speed_func_3d, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 1.0);
  }
  {
    olim6_rhr_arma m {n, 1, n, h, default_speed_func_3d, 0, 1, 1};
    m.add_boundary_node(1, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 1), 0.0);
  }

  /**
   * Tests for quadrants in the y-z plane:
   */
  {
    olim6_rhr_arma m {1, n, n, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0 + sqrt(2)/2.0);
  }
  {
    olim6_rhr_arma m {1, n, n, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim6_rhr_arma m {1, n, n, h, default_speed_func_3d, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim6_rhr_arma m {1, n, n, h, default_speed_func_3d, 1, 0, 1};
    m.add_boundary_node(0, 1, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0 + sqrt(2)/2.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 0.0);
  }
}

void octants_are_correct() {
  int n = 2;
  double h = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      for (int k = 0; k < 2; ++k) {
        double x0 = j, y0 = i, z0 = k;
        olim6_rhr_arma m {n, n, n, h, default_speed_func_3d, x0, y0, z0};
        m.add_boundary_node(i, j, k);
        m.run();
        IS_APPROX_EQUAL(m.get_value(i, j, k), 0.0);
        IS_APPROX_EQUAL(m.get_value((i + 1) % 2, j, k), 1.0);
        IS_APPROX_EQUAL(m.get_value(i, (j + 1) % 2, k), 1.0);
        IS_APPROX_EQUAL(m.get_value(i, j, (k + 1) % 2), 1.0);
        IS_APPROX_EQUAL(
          m.get_value((i + 1) % 2, (j + 1) % 2, k), 1.0 + sqrt(2)/2.0);
        IS_APPROX_EQUAL(
          m.get_value((i + 1) % 2, j, (k + 1) % 2), 1.0 + sqrt(2)/2.0);
        IS_APPROX_EQUAL(
          m.get_value(i, (j + 1) % 2, (k + 1) % 2), 1.0 + sqrt(2)/2.0);
        IS_APPROX_EQUAL(
          m.get_value((i + 1) % 2, (j + 1) % 2, (k + 1) % 2),
          1.0 + sqrt(2)/2.0 + sqrt(3)/3.0);
      }
    }
  }
}

void planes_are_correct() {
  int n = 11;
  double h = 1.0/(n/2);
  
  olim4_rhr m4 {n, n, h, default_speed_func, 1, 1};
  m4.add_boundary_node(n/2, n/2);
  m4.run();
  
  olim6_rhr_arma m6 {n, n, n, h, default_speed_func_3d, 1, 1, 1};
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
  olim6_rhr_arma m {n, n, n, 1.0, default_speed_func_3d, 1.0, 1.0, 1.0};
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
  quadrants_are_correct();
  octants_are_correct();
  planes_are_correct();
  result_is_symmetric();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
