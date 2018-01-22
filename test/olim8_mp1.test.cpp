#include "olim.test.common.hpp"

#include <olim.hpp>
#include <src/config.hpp>

int main() {
  using olim = olim8_mp1;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
  quadrants_are_correct<olim>(sqrt(2));
  {
#if TRIAL_NODE_OPTIMIZATION
    double tol = EPS(double);
#else
    double tol = 10*EPS(double);
#endif
    correct_corners_in_limit<olim>(101, tol);
  }
  result_is_symmetric<olim>((speed_func) default_speed_func);
  result_is_symmetric<olim>((speed_func) s1);
}
