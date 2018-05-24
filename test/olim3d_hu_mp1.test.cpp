#include "olim3d.hpp"
#include "olim.test.common.hpp"

using olim3d_t = olim3d_hu_mp1;

TEST (olim3d_hu_mp1, quadrants_are_correct) {
  quadrants_are_correct<olim3d_t>(sqrt(2));
}

