#include "olim.hpp"
#include "olim3d.hpp"
#include "olim.test.common.hpp"

#define MARCHER olim6_mp0

using olim_t = olim4_mp0;
using olim3d_t = MARCHER;

TEST (MARCHER, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(1.0 + sqrt(2)/2));
}

TEST (MARCHER, octants_are_correct) {
  ASSERT_TRUE(octants_are_correct<olim3d_t>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3));
}

TEST (MARCHER, planes_are_correct) {
  auto res = planes_are_correct<olim_t, olim3d_t>(speed_funcs[0], speed_funcs_3d[0]);
  ASSERT_TRUE(res);

  res = planes_are_correct<olim_t, olim3d_t>(speed_funcs[1], speed_funcs_3d[1], 9);
  ASSERT_TRUE(res);
}

TEST (MARCHER, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim3d_t>(speed_funcs_3d[0]));
  ASSERT_TRUE(result_is_symmetric<olim3d_t>(speed_funcs_3d[1]));
}

TEST (MARCHER, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (MARCHER, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}

TEST (MARCHER, s1_symmetry_test) {
  int n = 21;
  double h = 2./(n - 1);
  int i0 = n/2;

  olim3d_t o {n, n, n, h, speed_funcs_3d[1], 1., 1., 1.};
  o.add_boundary_node(i0, i0, i0);
  o.run();

  int i = 5, j1 = 3, j2 = 17, k = 9;
  double u = 0.5590785434269516;
  double U1 = o.get_value(i, j1, k);
  double U2 = o.get_value(i, j2, k);

  ASSERT_NEAR(U1, u, 2e-16);
  ASSERT_NEAR(U2, u, 2e-16);
}

TEST (MARCHER, solution_is_exact_in_factored_region) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim_t>(3));
}
