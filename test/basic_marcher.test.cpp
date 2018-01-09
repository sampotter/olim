#include "basic_marcher.hpp"
#include "olim.test.common.hpp"

int main() {
  using olim = basic_marcher;

  trivial_case_works<olim>();
  adjacent_update_works<olim>();
  correct_corners_in_limit<olim>(101, 1.86e-2);
}

