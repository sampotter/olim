#include "olim.hpp"
#include "olim3d.hpp"
#include "olim.test.common.hpp"

using olim_t = olim8_mp0;
using olim3d_t = olim26_mp0;

TEST (olim26_mp0, quadrants_are_correct) {
  quadrants_are_correct<olim3d_t>(sqrt(2));
}

TEST (olim26_mp0, octants_are_correct) {
  octants_are_correct<olim3d_t>(sqrt(2), sqrt(3));
}

TEST (olim26_mp0, planes_are_correct) {
  for (int i = 0; i < 2; ++i) {
    planes_are_correct<olim_t, olim3d_t>(speed_funcs[i], speed_funcs_3d[i]);
  }
}

TEST (olim26_mp0, result_is_symmetric) {
  for (int i = 0; i < 2; ++i) {
    result_is_symmetric<olim3d_t>(speed_funcs_3d[i]);
  }
}

TEST (olim26_mp0, two_by_two_by_three_cells_are_correct) {
  two_by_two_by_three_cells_are_correct<olim3d_t>();
}

TEST (olim26_mp0, plane_boundaries_are_correct) {
  plane_boundaries_are_correct<olim3d_t>();
}
