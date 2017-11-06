#include "olim.test.common.hpp"
#include "olim26.hpp"
#include "olim8.hpp"

int main() {
  quadrants_are_correct<olim26_mp1>(sqrt(2));
  octants_are_correct<olim26_mp1>(sqrt(2), sqrt(3));
  planes_are_correct<olim8_mp1, olim26_mp1>();
}
