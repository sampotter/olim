#include "test.hpp"

#include "basic_marcher_3d.hpp"
#include "olim.test.common.hpp"
#include "olim4.hpp"
#include "olim6.hpp"

int main() {
  quadrants_are_correct<olim6_mp1>(1 + sqrt(2)/2);
  octants_are_correct<olim6_mp1>(1.0 + sqrt(2)/2, 1.0 + sqrt(2)/2 + sqrt(3)/3);
  planes_are_correct<olim4_mp1, olim6_mp1>();
}
