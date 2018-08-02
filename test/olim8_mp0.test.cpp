#include "olim.hpp"
#include "olim.test.common.hpp"

using olim_t = olim8_mp0;

TEST (olim8_mp0, trivial_case_works) {
  ASSERT_TRUE(trivial_case_works<olim_t>());
}

TEST (olim8_mp0, adjacent_update_works) {
  ASSERT_TRUE(adjacent_update_works<olim_t>());
}

TEST (olim8_mp0, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim_t>(sqrt(2)));
}

TEST (olim8_mp0, correct_corners_in_limit) {
  ASSERT_TRUE(correct_corners_in_limit<olim_t>(101, 10*EPS(double)));
}

TEST (olim8_mp0, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) default_speed_func));
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) s1));
}
