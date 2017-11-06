#ifndef __OLIM6_IMPL_HPP__
#define __OLIM6_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"
#include "olim.macros.hpp"

template <class line_updates, class tri_updates, class tetra_updates>
void olim6<line_updates, tri_updates, tetra_updates>::update_impl(
  int i, int j, int k, double & T)
{
  using std::min;

#ifdef PRINT_UPDATES
  printf("olim6::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  abstract_node * nb[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  get_valid_neighbors(i, j, k, nb);
  double h = get_h(), s = speed(i, j, k), s_[6];
  for (int l = 0; l < 6; ++l) {
    if (nb[l]) {
      s_[l] = speed(i + di[l], j + dj[l], k + dk[l]);
    }
  }

  for (int l0 = 0, l1 = 1, l2 = 2;
       l0 < 6;
       ++l0, l1 = (l1 + 1) % 6, l2 = (l2 + 1) % 6) {
    if (nb[l0]) {
      LINE1(l0);
      if (nb[l1]) TRI11(l0, l1);
      if (nb[l2]) TRI11(l0, l2);
      if (nb[l1] && nb[l2]) TETRA111(l0, l1, l2);
    }
  }
  if (nb[0] && nb[2] && nb[4]) TETRA111(0, 2, 4);
  if (nb[1] && nb[3] && nb[5]) TETRA111(1, 3, 5);

#ifdef PRINT_UPDATES
  printf("olim6::update_impl: T <- %g\n", T);
#endif
}

#endif // __OLIM6_IMPL_HPP__
