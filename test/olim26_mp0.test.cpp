#include "olim.hpp"
#include "olim3d.hpp"
#include "olim.test.common.hpp"

#define MARCHER olim26_mp0

using olim_t = olim8_mp0;
using olim3d_t = MARCHER;

TEST (MARCHER, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim3d_t>(sqrt(2)));
}

TEST (MARCHER, octants_are_correct) {
  ASSERT_TRUE(octants_are_correct<olim3d_t>(sqrt(2), sqrt(3)));
}

TEST (MARCHER, two_by_two_by_three_cells_are_correct) {
  ASSERT_TRUE(two_by_two_by_three_cells_are_correct<olim3d_t>());
}

TEST (MARCHER, plane_boundaries_are_correct) {
  ASSERT_TRUE(plane_boundaries_are_correct<olim3d_t>());
}

TEST (MARCHER, solution_is_exact_in_factored_region) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim3d_t>(5));
}
