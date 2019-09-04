#include "common.hpp"

#include "eikonal/olim.hpp"
#include "eikonal/olim3d.hpp"

using namespace eikonal;

using olim_t = olim8_mp1<>;
using olim3d_t = olim18_mp1<>;

TEST (olim18_mp1, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(sqrt(2)));
}

TEST (olim18_mp1, octants_are_correct) {
  ASSERT_TRUE(octants_are_correct<olim3d_t>(sqrt(2), sqrt(2) + 1.0/sqrt(3)));
}

TEST (olim18_mp1, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (olim18_mp1, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}

TEST (olim18_mp1, solution_is_exact_in_factored_square) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim3d_t>(5));
}
