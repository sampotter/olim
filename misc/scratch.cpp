#include <olim3d.hpp>

int main() {
  int n = 11;
  olim18_rhr m {n, n, n, 2./(n-1), (speed_func_3d) default_speed_func, 1., 1., 1.};
  m.add_boundary_node(n/2, n/2, n/2);
  m.run();
}
