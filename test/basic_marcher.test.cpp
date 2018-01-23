#include "basic_marcher.hpp"
#include "olim.test.common.hpp"

using olim = basic_marcher;

TEST (basic_marcher, trivial_case_works) {
  trivial_case_works<olim>();
}

TEST (basic_marcher, adjacent_update_works) {
  adjacent_update_works<olim>();
}

TEST (basic_marcher, quadrants_are_correct) {
  quadrants_are_correct<olim>(1.0 + 1.0/sqrt(2));
}

TEST (basic_marcher, correct_corners_in_limit) {
  correct_corners_in_limit<olim>(101, 2.7e-2);
}

TEST (basic_marcher, result_is_symmetric) {
  result_is_symmetric<olim>((speed_func) default_speed_func);
  result_is_symmetric<olim>((speed_func) s1);
}

