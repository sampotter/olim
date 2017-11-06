#include "olim.test.common.hpp"
#include "olim26.hpp"
#include "olim8.hpp"

int main() {
  using olim = olim8_mp0;
  using olim3d = olim26_mp0;

  quadrants_are_correct<olim3d>(sqrt(2));
  octants_are_correct<olim3d>(sqrt(2), sqrt(3));
  planes_are_correct<olim, olim3d>();
}
