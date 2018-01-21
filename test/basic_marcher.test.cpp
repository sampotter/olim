#include "basic_marcher.hpp"
#include "olim.test.common.hpp"

int main() {
  using olim = basic_marcher;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
  quadrants_are_correct<olim>(1.0 + 1.0/sqrt(2));
  correct_corners_in_limit<olim>(101, 1.86e-2);
  result_is_symmetric<olim>((speed_func) default_speed_func);
  result_is_symmetric<olim>((speed_func) s1);
}

