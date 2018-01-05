#include "basic_marcher.hpp"
#include "basic_marcher_3d.hpp"
#include "olim.test.common.hpp"
#include "test.hpp"

int main() {
  using olim = basic_marcher;
  using olim3d = basic_marcher_3d;

  quadrants_are_correct<olim3d>(1.0 + sqrt(2)/2);
  octants_are_correct<olim3d>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3);
  planes_are_correct<olim, olim3d>();
  result_is_symmetric<olim3d>((speed_func_3d) default_speed_func);
  result_is_symmetric<olim3d>((speed_func_3d) s1);
  two_by_two_by_three_cells_are_correct<olim3d>();
  plane_boundaries_are_correct<olim3d>();
}
