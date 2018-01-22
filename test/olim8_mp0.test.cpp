#include "olim.test.common.hpp"

#include <olim.hpp>
#include <src/config.hpp>

int main() {
  using olim = olim8_mp0;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
  quadrants_are_correct<olim>(sqrt(2));
  correct_corners_in_limit<olim>(101, 10*EPS(double));
  result_is_symmetric<olim>((speed_func) default_speed_func);
  result_is_symmetric<olim>((speed_func) s1);
}
