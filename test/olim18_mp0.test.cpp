#include "olim.hpp"
#include "olim3d.hpp"
#include "olim.test.common.hpp"

using olim_t = olim8_mp0;
using olim3d_t = olim18_mp0;

TEST (olim18_mp0, quadrants_are_correct) {
  quadrants_are_correct<olim3d_t>(sqrt(2));
}

TEST (olim18_mp0, octants_are_correct) {
  octants_are_correct<olim3d_t>(sqrt(2), sqrt(2) + 1.0/sqrt(3));
}

TEST (olim18_mp0, planes_are_correct) {
  planes_are_correct<olim_t, olim3d_t>(speed_funcs[0], speed_funcs_3d[0]);
  planes_are_correct<olim_t, olim3d_t>(speed_funcs[1], speed_funcs_3d[1]);
}

TEST (olim18_mp0, result_is_symmetric) {
  result_is_symmetric<olim3d_t>(speed_funcs_3d[0]);
  result_is_symmetric<olim3d_t>(speed_funcs_3d[1]);
}

TEST (olim18_mp0, two_by_two_by_three_cells_are_correct) {
  two_by_two_by_three_cells_are_correct<olim3d_t>();
}

TEST (olim18_mp0, plane_boundaries_are_correct) {
  plane_boundaries_are_correct<olim3d_t>();
}

TEST (olim18_mp0, planes_are_correct_nonsymmetric) {
  planes_are_correct_nonsymmetric<olim_t, olim3d_t>(s4, f4xy, f4yz, f4xz);
}

TEST (olim18_mp0, planes_agree_nonsymmetric) {
  planes_agree_nonsymmetric<olim_t, olim3d_t>(s4, s4xy, s4yz, s4xz);
}
