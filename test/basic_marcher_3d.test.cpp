#include "basic_marcher.hpp"
#include "basic_marcher_3d.hpp"
#include "olim.test.common.hpp"
#include "test.hpp"

void plane_boundaries_are_correct() {
  int n = 2;
  double h = 1;
  basic_marcher_3d m {n, n, n, h, default_speed_func_3d, 0, 0, 0};

  node_3d nodes[4];
  for (int i = 0, k = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      nodes[k++] = node_3d {i, j, 0};
    }
  }
  m.add_boundary_nodes(nodes, 4);

  m.run();
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 2; ++j) {
      IS_APPROX_EQUAL(m.get_value(i, j, 1), 1.0);
    }
  }
}

int main() {
  using olim = basic_marcher;
  using olim3d = basic_marcher_3d;

  quadrants_are_correct<olim3d>(1.0 + sqrt(2)/2);
  octants_are_correct<olim3d>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3);
  planes_are_correct<olim, olim3d>();
  result_is_symmetric<olim3d>();
  two_by_two_by_three_cells_are_correct<olim3d>();
  plane_boundaries_are_correct();
}
