#include "test.hpp"

#include "basic_marcher_3d.hpp"
#include "olim.test.common.hpp"
#include "olim4.hpp"
#include "olim6.hpp"

int main() {
  using olim = olim4_mp1;
  using olim3d = olim6_mp1;

  quadrants_are_correct<olim3d>(1 + sqrt(2)/2);
  octants_are_correct<olim3d>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3);
  planes_are_correct<olim, olim3d>();
  two_by_two_by_three_cells_are_correct<olim3d>();
  result_is_symmetric<olim3d>();
  plane_boundaries_are_correct<olim3d>();

#if RELWITHDEBINFO
  int n = 31;
#else
  int n = 11;
#endif
  agrees_with_other_olim3d<olim3d, olim6_rhr>(n);
  agrees_with_other_olim3d<olim3d, olim6_mp0>(n);
}
