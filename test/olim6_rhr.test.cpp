#include "olim.test.common.hpp"

#include <basic_marcher_3d.hpp>
#include <olim.hpp>
#include <olim3d.hpp>

void agrees_with_basic_marcher_3d(speed_func_3d s) {
  int n = 11;
  double h = 2.0/(n - 1);
  int i0 = (n - 1)/2, j0 = i0, k0 = i0;

  basic_marcher_3d m3d {n, n, n, h, s, 1, 1, 1};
  m3d.add_boundary_node(i0, j0, k0);
  m3d.run();

  olim6_rhr m6 {n, n, n, h, s, 1, 1, 1};
  m6.add_boundary_node(i0, j0, k0);
  m6.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        IS_APPROX_EQUAL(m3d.get_value(i, j, k), m6.get_value(i, j, k), 1e-13);
      }
    }
  }
}

int main() {
  using olim = olim4_rhr;
  using olim3d = olim6_rhr;

  quadrants_are_correct<olim3d>(1 + sqrt(2)/2);
  octants_are_correct<olim3d>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3);
  for (int i = 0; i < 2; ++i) {
    planes_are_correct<olim, olim3d>(speed_funcs[i], speed_funcs_3d[i]);
    result_is_symmetric<olim3d>(speed_funcs_3d[i]);
    agrees_with_basic_marcher_3d(speed_funcs_3d[i]);
  }
  two_by_two_by_three_cells_are_correct<olim3d>();
  plane_boundaries_are_correct<olim3d>();
}
