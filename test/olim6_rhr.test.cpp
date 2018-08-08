#include "olim.test.common.hpp"

#include <basic_marcher_3d.hpp>
#include <olim.hpp>
#include <olim3d.hpp>

#define MARCHER olim6_rhr

using olim_t = olim4_rhr;
using olim3d_t = olim6_rhr;

testing::AssertionResult
agrees_with_basic_marcher_3d(speed_func_3d s) {
  int n = 11;
  double h = 2.0/(n - 1);
  int i0 = (n - 1)/2, j0 = i0, k0 = i0;

  basic_marcher_3d m3d {n, n, n, h, s, 1, 1, 1};
  m3d.add_boundary_node(i0, j0, k0);
  m3d.run();

  olim3d_t m6 {n, n, n, h, s, 1, 1, 1};
  m6.add_boundary_node(i0, j0, k0);
  m6.run();

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        auto res = eq_w_rel_tol(m3d.get_value(i, j, k), m6.get_value(i, j, k));
        if (!res) return res;
      }
    }
  }

  return testing::AssertionSuccess();
}

TEST (MARCHER, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(1 + sqrt(2)/2));
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
  for (int i = 0; i < 2; ++i) {
    ASSERT_TRUE(result_is_symmetric<olim3d_t>(speed_funcs_3d[i]));
  }
}

TEST (MARCHER, agrees_with_basic_marcher_3d) {
  for (int i = 0; i < 2; ++i) {
    ASSERT_TRUE(agrees_with_basic_marcher_3d(speed_funcs_3d[i]));
  }
}

TEST (MARCHER, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (MARCHER, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}

TEST (MARCHER, solution_is_exact_in_factored_region) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim_t>(3));
}
