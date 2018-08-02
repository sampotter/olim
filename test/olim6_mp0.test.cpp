#include "olim.hpp"
#include "olim3d.hpp"
#include "olim.test.common.hpp"

using olim_t = olim4_mp0;
using olim3d_t = olim6_mp0;

TEST (olim6_mp0, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(1.0 + sqrt(2)/2));
}

TEST (olim6_mp0, octants_are_correct) {
  ASSERT_TRUE(octants_are_correct<olim3d_t>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3));
}

TEST (olim6_mp0, planes_are_correct) {
  auto res = planes_are_correct<olim_t, olim3d_t>(speed_funcs[0], speed_funcs_3d[0]);
  ASSERT_TRUE(res);

  res = planes_are_correct<olim_t, olim3d_t>(speed_funcs[1], speed_funcs_3d[1], 9);
  ASSERT_TRUE(res);
}

TEST (olim6_mp0, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim3d_t>(speed_funcs_3d[0]));
  ASSERT_TRUE(result_is_symmetric<olim3d_t>(speed_funcs_3d[1]));
}

TEST (olim6_mp0, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (olim6_mp0, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}

TEST (olim6_mp0, s1_symmetry_test) {
  int n = 21;
  double h = 2./(n - 1);
  int i0 = n/2;

  olim3d_t o {n, n, n, h, speed_funcs_3d[1], 1., 1., 1.};
  o.add_boundary_node(i0, i0, i0);
  o.run();

  int i = 5, j1 = 3, j2 = 17, k = 9;
  double u = 0.5590785434269516;

  ASSERT_EQ(o.get_value(i, j1, k), u);
  ASSERT_EQ(o.get_value(i, j2, k), u);
}
