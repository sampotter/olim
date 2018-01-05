#include "olim.test.common.hpp"

#include <olim.hpp>
#include <olim3d.hpp>

int main() {
  using olim = olim8_mp0;
  using olim3d = olim26_mp0;

  quadrants_are_correct<olim3d>(sqrt(2));
  octants_are_correct<olim3d>(sqrt(2), sqrt(2) + 1.0/sqrt(3));
  for (int i = 0; i < 2; ++i) {
    planes_are_correct<olim, olim3d>(speed_funcs[i], speed_funcs_3d[i]);
    result_is_symmetric<olim3d>(speed_funcs_3d[i]);
  }
  two_by_two_by_three_cells_are_correct<olim3d>();
  plane_boundaries_are_correct<olim3d>();
}
