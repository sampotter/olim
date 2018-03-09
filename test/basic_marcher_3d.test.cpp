#include "basic_marcher.hpp"
#include "basic_marcher_3d.hpp"
#include "olim.test.common.hpp"

using olim = basic_marcher;
using olim3d = basic_marcher_3d;

TEST (basic_marcher_3d, quadrants_are_correct) {
  quadrants_are_correct<olim3d>(1.0 + sqrt(2)/2);
}

TEST (basic_marcher_3d, octants_are_correct) {
  octants_are_correct<olim3d>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3);
}

TEST (basic_marcher_3d, planes_are_correct) {
  planes_are_correct<olim, olim3d>();
}

TEST (basic_marcher_3d, result_is_symmetric) {
  result_is_symmetric<olim3d>((speed_func_3d) default_speed_func);
  result_is_symmetric<olim3d>((speed_func_3d) s1);
}

TEST (basic_marcher_3d, two_by_two_by_three_cells_are_correct) {
  two_by_two_by_three_cells_are_correct<olim3d>();
}

TEST (basic_marcher_3d, plane_boundaries_are_correct) {
  plane_boundaries_are_correct<olim3d>();
}

TEST (basic_marcher_3d, planes_are_correct_nonsymmetric) {
  planes_are_correct_nonsymmetric<olim, olim3d>(s4, s4xy, s4yz, s4xz);
}

TEST (basic_marcher_3d, planes_agree_nonsymmetric) {
  planes_agree_nonsymmetric<olim, olim3d>(s4, s4xy, s4yz, s4xz);
}
