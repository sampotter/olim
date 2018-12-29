#include "olim.test.common.hpp"
#include "olim.hpp"

using olim_t = olim8_mp1;

TEST (olim8_mp1, trivial_case_works) {
  ASSERT_TRUE(trivial_case_works<olim_t>());
}

TEST (olim8_mp1, adjacent_update_works) {
  ASSERT_TRUE(adjacent_update_works<olim_t>());
}

TEST (olim8_mp1, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim_t>(sqrt(2)));
}

TEST (olim8_mp1, correct_corners_in_limit) {
  ASSERT_TRUE(correct_corners_in_limit<olim_t>(101, 10*eps<double>));
}

TEST (olim8_mp1, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) default_speed_func));
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) s1));
}

TEST (olim8_mp1, solution_is_exact_in_factored_region) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim_t>(5));
}
