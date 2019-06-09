#include "common.hpp"

#include "eikonal/fmm.hpp"
#include "eikonal/olim.hpp"
#include "eikonal/olim3d.hpp"

using olim_t = olim4_rhr<>;
using olim3d_t = olim6_rhr<>;

TEST (olim6_rhr, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(1 + sqrt(2)/2));
}

TEST (olim6_rhr, octants_are_correct) {
  ASSERT_TRUE(octants_are_correct<olim3d_t>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3));
}

TEST (olim6_rhr, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (olim6_rhr, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}

TEST (olim6_rhr, solution_is_exact_in_factored_region) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim3d_t>(5));
}
