#include "olim.test.common.hpp"
#include "olim.hpp"

using olim_t = olim8_rhr;

TEST (olim8_rhr, trivial_case_works) {
  trivial_case_works<olim_t>();
}

TEST (olim8_rhr, adjacent_update_works) {
  adjacent_update_works<olim_t>();
}

TEST (olim8_rhr, quadrants_are_correct) {
  quadrants_are_correct<olim_t>(sqrt(2));
}

TEST (olim8_rhr, correct_corners_in_limit) {
  correct_corners_in_limit<olim_t>(101, 10*EPS(double));
}

TEST (olim8_rhr, result_is_symmetric) {
  result_is_symmetric<olim_t>((speed_func) default_speed_func);
  result_is_symmetric<olim_t>((speed_func) s1);
}
