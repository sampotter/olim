#include "basic_marcher.hpp"
#include "basic_marcher_3d.hpp"
#include "olim.test.common.hpp"

using olim = basic_marcher;
using olim3d_t = basic_marcher_3d;

TEST (basic_marcher_3d, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(1.0 + sqrt(2)/2));
}

TEST (basic_marcher_3d, octants_are_correct) {
  ASSERT_TRUE(
    octants_are_correct<olim3d_t>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3));
}

TEST (basic_marcher_3d, planes_are_correct) {
  for (auto i = 0ull; i < std::size(speed_funcs_3d); ++i) {
    auto s = speed_funcs[i];
    auto s3d = speed_funcs_3d[i];
    auto res = planes_are_correct<olim, olim3d_t>(s, s3d);
    ASSERT_TRUE(res);
  }
}

TEST (basic_marcher_3d, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim3d_t>((speed_func_3d) default_speed_func));
  ASSERT_TRUE(result_is_symmetric<olim3d_t>((speed_func_3d) s1, 21, 1e-10));
}

TEST (basic_marcher_3d, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (basic_marcher_3d, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}
