#include "olim.test.common.hpp"

#include <olim.hpp>

using olim_t = olim4_mp1;

TEST (olim4_mp1, trivial_case_works) {
  trivial_case_works<olim_t>();
}

TEST (olim4_mp1, adjacent_update_works) {
  adjacent_update_works<olim_t>();
}

TEST (olim4_mp1, quadrants_are_correct) {
  quadrants_are_correct<olim_t>(1.0 + 1.0/sqrt(2));
}

TEST (olim4_mp1, correct_corners_in_limit) {
  correct_corners_in_limit<olim_t>(101, 2.7e-2);
}

TEST (olim4_mp1, result_is_symmetric) {
  result_is_symmetric<olim_t>((speed_func) default_speed_func);
  result_is_symmetric<olim_t>((speed_func) s1);
}
