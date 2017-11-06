#ifndef __OLIM4_IMPL_HPP__
#define __OLIM4_IMPL_HPP__

#include <algorithm>

#include "common.macros.hpp"
#include "olim.macros.hpp"

template <class line_updates, class tri_updates>
void olim4<line_updates, tri_updates>::update_impl(int i, int j, double & T) {
  using std::min;

  abstract_node * nb[4] = {0x0, 0x0, 0x0, 0x0}; // NESW
  get_valid_neighbors(i, j, nb);
  double s = speed(i, j), h = get_h(), s_[4];
  for (int k = 0; k < 4; ++k) {
    if (nb[k]) {
      s_[k] = speed(i + di[k], j + dj[k]);
    }
  }

  for (int k = 0; k < 4; ++k) {
    if (nb[k] && !nb[(k + 1) % 4] && !nb[(k + 3) % 4]) {
      LINE1(k);
    }
  }

  for (int k = 0, l = k + 1; k < 4; ++k, l = (k + 1) % 4) {
    if (nb[k] && nb[l]) {
      TRI11(k, l);
    }
  }
}

#endif // __OLIM4_IMPL_HPP__
