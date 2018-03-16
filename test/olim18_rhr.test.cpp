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
  planes_are_correct<olim_t, olim3d_t>(speed_funcs[0], speed_funcs_3d[0]);

  // TODO: again---some kind of CFL condition might be being violated
  // for n < 9... would depend on method and speed function
  planes_are_correct<olim_t, olim3d_t>(speed_funcs[1], speed_funcs_3d[1], 9);
}

TEST (olim18_rhr, result_is_symmetric) {
  result_is_symmetric<olim3d_t>(speed_funcs_3d[0]);
  result_is_symmetric<olim3d_t>(speed_funcs_3d[1]);
}

TEST (olim18_rhr, two_by_two_by_three_cells_are_correct) {
  two_by_two_by_three_cells_are_correct<olim3d_t>();
}

TEST (olim18_rhr, plane_boundaries_are_correct) {
  plane_boundaries_are_correct<olim3d_t>();
}

TEST (olim18_rhr, planes_are_correct_nonsymmetric) {
  planes_are_correct_nonsymmetric<olim_t, olim3d_t>(s4, f4xy, f4yz, f4xz);
}

TEST (olim18_rhr, planes_agree_nonsymmetric) {
  planes_agree_nonsymmetric<olim_t, olim3d_t>(s4, s4xy, s4yz, s4xz);
}
