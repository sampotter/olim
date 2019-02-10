#include "fmm.hpp"
#include "common.hpp"
#include "olim.hpp"
#include "olim.test.common.hpp"

using olim_t = olim4_rhr;

TEST (olim4_rhr, trivial_case_works) {
  ASSERT_TRUE(trivial_case_works<olim_t>());
}

TEST (olim4_rhr, adjacent_update_works) {
  ASSERT_TRUE(adjacent_update_works<olim_t>());
}

TEST (olim4_rhr, quadrants_are_correct) {
  ASSERT_TRUE(quadrants_are_correct<olim_t>(1.0 + 1.0/sqrt(2)));
}

TEST (olim4_rhr, correct_corners_in_limit) {
  ASSERT_TRUE(correct_corners_in_limit<olim_t>(101, 2.7e-2));
}

TEST (olim4_rhr, result_is_symmetric) {
  ASSERT_TRUE(result_is_symmetric<olim_t>(s0<2>));
  ASSERT_TRUE(result_is_symmetric<olim_t>(s1<2>));
}

TEST (olim4_rhr, agrees_with_basic_marcher) {
  for (auto & s: slow2s) {
    auto res = olims_agree<olim_t, fmm<2>>(s);
    ASSERT_TRUE(res);
  }
}

TEST (olim4_rhr, solution_is_exact_in_factored_region) {
  ASSERT_TRUE(solution_is_exact_in_factored_square<olim_t>(3));
}
