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

  // Check that planes are correct:
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double m8value = m8.get_value(i, j);
      IS_APPROX_EQUAL(m8value, m26.get_value(i, j, n/2));
      IS_APPROX_EQUAL(m8value, m26.get_value(i, n/2, j));
      IS_APPROX_EQUAL(m8value, m26.get_value(n/2, i, j));
    }
  }
}

int main() {
  planes_are_correct();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
