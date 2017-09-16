#include "test.hpp"

#include "olim4_rhr.hpp"
#include "olim6_rhr.hpp"

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

int main() {
  planes_are_correct();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
