#include "basic_marcher.hpp"
#include "test_graph_marcher.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>

int main() {
  int n = 21;
  double h = 2.0/(n - 1);

  test_graph_marcher m1 {n, n, h, default_speed_func, 1, 1};
  m1.add_boundary_node((n - 1)/2, (n - 1)/2);
  m1.run();

  basic_marcher m2 {n, n, h, default_speed_func, 1, 1};
  m2.add_boundary_node((n - 1)/2, (n - 1)/2);
  m2.run();

  double maxerr {0};
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      double u1 = m1.get_value(i, j);
      double u2 = m2.get_value(i, j);
      maxerr = std::max(maxerr, std::fabs(u1 - u2));
    }
  }

  printf("maxerr = %g\n", maxerr);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
