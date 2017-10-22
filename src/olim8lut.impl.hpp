#ifndef __OLIM8LUT_IMPL_HPP__
#define __OLIM8LUT_IMPL_HPP__

#include <algorithm>

#include <src/config.hpp>

#include "common.macros.hpp"
#include "olim.macros.def.hpp"

template <class line_updates, class tri_updates>
void olim8lut<line_updates, tri_updates>::update_impl(int i, int j, double & T) {
  using std::min;

  abstract_node * nb[8] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  get_valid_neighbors(i, j, nb);
  double s = speed(i, j), h = get_h();

  // implementing constrained algorithm for now
  for (int k = 0, l = 1, m = 2; k < 8; k += 2, l = k + 1, m = (l + 1) % 8) {
    switch (!!nb[k] + 2*!!nb[l] + 2*!!nb[m]) {
    case 2:
      LINE2(k);
      break;
    case 4:
      TRI12(l, m);
      break;
    case 3:
    case 6:
      LINE1(k);
      break;
    case 5: 
    case 7:
      TRI12(k, l);
      break;
    default:
      break;
    }
  }
}

#include "olim.macros.undef.hpp"

#endif // __OLIM8LUT_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
