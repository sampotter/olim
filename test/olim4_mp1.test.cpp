#include "olim.test.common.hpp"

#include <olim.hpp>

using olim_t = olim4_mp1;

TEST (olim4_mp1, trivial_case_works) {
  ASSERT_TRUE(trivial_case_works<olim_t>());
}

TEST (olim4_mp1, adjacent_update_works) {
  ASSERT_TRUE(adjacent_update_works<olim_t>());
}

TEST (olim4_mp1, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim_t>(1.0 + 1.0/sqrt(2)));
}

TEST (olim4_mp1, correct_corners_in_limit) {
  ASSERT_TRUE(correct_corners_in_limit<olim_t>(101, 2.7e-2));
}

TEST (olim4_mp1, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim_t>(s0<2>));
  ASSERT_TRUE(result_is_symmetric<olim_t>(s1<2>));
}

TEST (olim4_mp1, solution_is_exact_in_factored_region) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim_t>(3));
}
