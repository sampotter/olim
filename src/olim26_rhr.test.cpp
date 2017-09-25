#include "test.hpp"

#include <cmath>

#include "olim8_rhr.hpp"
#include "olim26_rhr.hpp"

void planes_are_correct() {
  int n = 21;
  double h = 1.0/(n/2);
  
  olim8_rhr m8 {n, n, h, default_speed_func, 1, 1};
  m8.add_boundary_node(n/2, n/2);
  m8.run();
  
  olim26_rhr_arma m26 {n, n, n, h, default_speed_func_3d, 1, 1, 1};
  m26.add_boundary_node(n/2, n/2, n/2);
  m26.run();

  /**
   * Tests for quadrants in the y-z plane:
   */
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 0, 0, 0};
    m.add_boundary_node(0, 0, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), sqrt(2));
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 1, 0, 0};
    m.add_boundary_node(0, 1, 0);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 0, 0, 1};
    m.add_boundary_node(0, 0, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 0.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 1.0);
  }
  {
    olim26_rhr_arma m {1, n, n, h, default_speed_func_3d, 1, 0, 1};
    m.add_boundary_node(0, 1, 1);
    m.run();
    IS_APPROX_EQUAL(m.get_value(0, 0, 0), sqrt(2));
    IS_APPROX_EQUAL(m.get_value(0, 0, 1), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 0), 1.0);
    IS_APPROX_EQUAL(m.get_value(0, 1, 1), 0.0);
  }
}

int main() {
  planes_are_correct();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
