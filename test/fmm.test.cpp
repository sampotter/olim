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


TEST (fmm3, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(1.0 + sqrt(2)/2));
}

TEST (fmm3, octants_are_correct) {
  ASSERT_TRUE(
    octants_are_correct<olim3d_t>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3));
}

TEST (fmm3, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (fmm3, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}
