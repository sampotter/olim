#include "olim.test.common.hpp"
#include "olim8.hpp"
#include "olim18.hpp"

int main() {
  using olim = olim8_mp0;
  using olim3d = olim18_mp0;

  quadrants_are_correct<olim3d>(sqrt(2));
  octants_are_correct<olim3d>(sqrt(2), sqrt(2) + 1.0/sqrt(3));
  planes_are_correct<olim, olim3d>();
  result_is_symmetric<olim3d>();
  two_by_two_by_three_cells_are_correct<olim3d>();
  plane_boundaries_are_correct<olim3d>();
}
