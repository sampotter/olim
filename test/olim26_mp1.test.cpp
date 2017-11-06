#include "olim.test.common.hpp"
#include "olim26.hpp"
#include "olim8.hpp"

int main() {
  using olim = olim8_mp1;
  using olim3d = olim26_mp1;

  quadrants_are_correct<olim3d>(sqrt(2));
  octants_are_correct<olim3d>(sqrt(2), sqrt(3));
  planes_are_correct<olim, olim3d>();
  result_is_symmetric<olim3d>();
  two_by_two_by_three_cells_are_correct<olim3d>();
}
