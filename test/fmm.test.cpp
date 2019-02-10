#include "fmm.hpp"
#include "olim.test.common.hpp"

using olim = fmm<2>;
using olim3d_t = fmm<3>;

TEST (fmm2, trivial_case_works) {
  ASSERT_TRUE(trivial_case_works<olim>());
}

TEST (fmm2, adjacent_update_works) {
  ASSERT_TRUE(adjacent_update_works<olim>());
}

TEST (fmm2, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim>(1.0 + 1.0/sqrt(2)));
}

TEST (fmm2, correct_corners_in_limit) {
  ASSERT_TRUE(correct_corners_in_limit<olim>(101, 2.7e-2));
}

TEST (fmm2, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim>(s0<2>));
  ASSERT_TRUE(result_is_symmetric<olim>(s1<2>));
}


TEST (fmm3, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(1.0 + sqrt(2)/2));
}

TEST (fmm3, octants_are_correct) {
  ASSERT_TRUE(
    octants_are_correct<olim3d_t>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3));
}

TEST (fmm3, planes_are_correct) {
  for (auto i = 0ull; i < std::size(slow3s); ++i) {
    auto s = slow2s[i];
    auto s3d = slow3s[i];
    auto res = planes_are_correct<olim, olim3d_t>(s, s3d);
    ASSERT_TRUE(res);
  }
}

TEST (fmm3, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim3d_t>(s0<3>));
  ASSERT_TRUE(result_is_symmetric<olim3d_t>(s1<3>, 21, 1e-10));
}

TEST (fmm3, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (fmm3, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}
