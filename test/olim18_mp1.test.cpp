#include "olim.test.common.hpp"
#include "olim18.hpp"
#include "olim8.hpp"

int main() {
  quadrants_are_correct<olim18_mp1>(sqrt(2));
  octants_are_correct<olim18_mp1>(sqrt(2), 1.0 + 2/sqrt(3));
  planes_are_correct<olim8_mp1, olim18_mp1>();
}
