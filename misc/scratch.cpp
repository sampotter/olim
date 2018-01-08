#include <olim.hpp>

int main() {
  int n = 1001;
  olim8_rhr m {n, n, 2./(n-1), (speed_func) default_speed_func, 1., 1.};
  m.add_boundary_node(50, 50);
  m.run();
}
