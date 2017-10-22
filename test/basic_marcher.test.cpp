#include "basic_marcher.hpp"
#include "speed_funcs.hpp"
#include "test.hpp"

void trivial_case_works() {
  basic_marcher m {1, 1};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 0.0);
}

void adjacent_update_works() {
  basic_marcher m {2, 1, 0.5};
  m.add_boundary_node(0, 0);
  m.run();
  IS_APPROX_EQUAL(m.get_value(1, 0), 0.5);
}

void neighboring_values_are_correct() {
  basic_marcher m {3, 3, 1};
  m.add_boundary_node(1, 1);
  m.run();
  double gt[] = {
    1.707106781186547,
    1,
    1.707106781186547,
    1,
    0,
    1,
    1.707106781186547,
    1,
    1.707106781186547
  };
  int k = 0;
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      IS_APPROX_EQUAL(gt[k++], m.get_value(i, j), 1e-15);
    }
  }
}

void rectangular_domain_works() {
  double h = 0.1, x0 = 1, y0 = 1;
  basic_marcher m {11, 21, h, default_speed_func, x0, y0};
  m.add_boundary_node(5, 10);
  m.run();
  IS_APPROX_EQUAL(m.get_value(0, 0), 1.17825, 1e-5);
}

int main() {
  trivial_case_works();
  adjacent_update_works();
  neighboring_values_are_correct();
  rectangular_domain_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
