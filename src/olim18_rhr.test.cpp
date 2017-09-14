#include "test.hpp"

#include "olim8_rhr.hpp"
#include "olim18_rhr.hpp"

void quadrants_are_correct() {
  int n = 2;
  double h = 1;
  {
    olim18_rhr_arma m {n, n, 1, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), sqrt(2));
  }
  {
    olim18_rhr_arma m {n, n, 1, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim18_rhr_arma m {n, n, 1, h, default_speed_func_3d, 0, 1, 0};
    m.add_boundary_node(1, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 1.0);
  }
  {
    olim18_rhr_arma m {n, n, 1, h, default_speed_func_3d, 1, 1, 0};
    m.add_boundary_node(1, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(1, 1, 0), 0.0);
  }
}

void planes_are_correct() {
  int n = 21;
  double h = 1.0/(n/2);
  
  olim8_rhr m8 {n, n, h, default_speed_func, 1, 1};
  m8.add_boundary_node(n/2, n/2);
  m8.run();
  
  olim18_rhr_arma m18 {n, n, n, h, default_speed_func_3d, 1, 1, 1};
  m18.add_boundary_node(n/2, n/2, n/2);
  m18.run();

  // Check that planes are correct:
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double m8value = m8.get_value(i, j);
      IS_APPROX_EQUAL(m8value, m18.get_value(i, j, n/2));
      IS_APPROX_EQUAL(m8value, m18.get_value(i, n/2, j));
      IS_APPROX_EQUAL(m8value, m18.get_value(n/2, i, j));
    }
  }
}

int main() {
  planes_are_correct();
  quadrants_are_correct();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
