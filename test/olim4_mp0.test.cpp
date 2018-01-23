#include "basic_marcher.hpp"
#include "olim.hpp"
#include "olim.test.common.hpp"

using olim_t = olim4_mp0;

TEST (olim4_mp0, trivial_case_works) {
  trivial_case_works<olim_t>();
}

TEST (olim4_mp0, adjacent_update_works) {
  adjacent_update_works<olim_t>();
}

TEST (olim4_mp0, quadrants_are_correct) {
  quadrants_are_correct<olim_t>(1.0 + 1.0/sqrt(2));
}

TEST (olim4_mp0, correct_corners_in_limit) {
  correct_corners_in_limit<olim_t>(101, 1.86e-2);
}

TEST (olim4_mp0, result_is_symmetric) {
  result_is_symmetric<olim_t>((speed_func) default_speed_func);
  result_is_symmetric<olim_t>((speed_func) s1);
}

TEST (olim4_rhr, agrees_with_basic_marcher) {
  olims_agree<olim_t, basic_marcher>((speed_func) default_speed_func);
}
