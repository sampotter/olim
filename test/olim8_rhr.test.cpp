#include "olim.test.common.hpp"
#include "olim.hpp"

using olim_t = olim8_rhr;

TEST (olim8_rhr, trivial_case_works) {
  ASSERT_TRUE(trivial_case_works<olim_t>());
}

TEST (olim8_rhr, adjacent_update_works) {
  ASSERT_TRUE(adjacent_update_works<olim_t>());
}

TEST (olim8_rhr, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim_t>(sqrt(2)));
}

TEST (olim8_rhr, correct_corners_in_limit) {
  ASSERT_TRUE(correct_corners_in_limit<olim_t>(101, 10*EPS(double)));
}

TEST (olim8_rhr, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) default_speed_func));
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) s1));
}

TEST (olim8_rhr, factoring_sanity_check) {
  ASSERT_TRUE(factoring_sanity_check<olim_t>((speed_func) s1, (speed_func) f1, 101));
}
