#include "common.defs.hpp"
#include "olim.hpp"
#include "olim3d.hpp"
#include "olim.test.common.hpp"

using olim_t = olim8_rhr;
using olim3d_t = olim18_rhr;

TEST (olim18_rhr, quadrants_are_correct) {
  quadrants_are_correct<olim3d_t>(sqrt2);
}

TEST (olim18_rhr, octants_are_correct) {
  octants_are_correct<olim3d_t>(sqrt2, sqrt2 + 1.0/sqrt3);
}

TEST (olim18_rhr, planes_are_correct) {
  for (int i = 0; i < 2; ++i) {
    planes_are_correct<olim_t, olim3d_t>(speed_funcs[i], speed_funcs_3d[i]);
  }
}

TEST (olim18_rhr, result_is_symmetric) {
  for (int i = 0; i < 2; ++i) {
    result_is_symmetric<olim3d_t>(speed_funcs_3d[i]);
  }
}

TEST (olim18_rhr, two_by_two_by_three_cells_are_correct) {
  two_by_two_by_three_cells_are_correct<olim3d_t>();
}

TEST (olim18_rhr, plane_boundaries_are_correct) {
  plane_boundaries_are_correct<olim3d_t>();
}
