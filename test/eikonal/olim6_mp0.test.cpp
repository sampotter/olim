#include "common.hpp"

#include "eikonal/olim.hpp"
#include "eikonal/olim3d.hpp"

using namespace eikonal;

using olim_t = olim4_mp0<>;
using olim3d_t = olim6_mp0<>;

TEST (olim6_mp0, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(1.0 + sqrt(2)/2));
}

TEST (olim6_mp0, octants_are_correct) {
  ASSERT_TRUE(octants_are_correct<olim3d_t>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3));
}

TEST (olim6_mp0, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (olim6_mp0, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}

TEST (olim6_mp0, solution_is_exact_in_factored_square) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim3d_t>(5));
}
