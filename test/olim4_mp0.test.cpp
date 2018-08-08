#include "basic_marcher.hpp"
#include "olim.hpp"
#include "olim.test.common.hpp"

using olim_t = olim4_mp0;

TEST (olim4_mp0, trivial_case_works) {
  ASSERT_TRUE(trivial_case_works<olim_t>());
}

TEST (olim4_mp0, adjacent_update_works) {
  ASSERT_TRUE(adjacent_update_works<olim_t>());
}

TEST (olim4_mp0, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim_t>(1.0 + 1.0/sqrt(2)));
}

TEST (olim4_mp0, correct_corners_in_limit) {
  ASSERT_TRUE(correct_corners_in_limit<olim_t>(101, 2.7e-2));
}

TEST (olim4_mp0, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) default_speed_func));
  ASSERT_TRUE(result_is_symmetric<olim_t>((speed_func) s1));
}

TEST (olim4_rhr, agrees_with_basic_marcher) {
  auto res = olims_agree<olim_t, basic_marcher>((speed_func) default_speed_func);
  ASSERT_TRUE(res);
}

TEST (olim4_mp0, solution_is_exact_in_factored_region) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim_t>(3));
}
