#ifndef __OLIM6_RECT_IMPL_HPP__
#define __OLIM6_RECT_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"
#include "olim6.defs.hpp"

template <class update_rules, class speed_estimate>
void olim6_rect<update_rules, speed_estimate>::update_impl(
  int i, int j, int k, double & T)
{
  using namespace olim6;
  using std::min;

#ifdef PRINT_UPDATES
  printf("olim6_rect::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  abstract_node * nb[6] = {0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
  get_valid_neighbors(i, j, k, nb);
  double h = get_h(), s = speed(i, j, k), s_hat, s_[6];
  for (int l = 0; l < 6; ++l) {
    if (nb[l]) {
      s_[l] = speed(i + di[l], j + dj[l], k + dk[l]);
    }
  }

  for (int l0 = 0, l1 = 1, l2 = 2; l0 < 6;
       ++l0, l1 = (l1 + 1) % 6, l2 = (l2 + 1) % 6) {
    if (nb[l0]) {
      s_hat = this->estimate_speed(s, s_[l0]);
      T = min(T, this->line1(VAL(l0), s_hat, h));
      if (nb[l1]) {
        s_hat = this->estimate_speed(s, s_[l0], s_[l1]);
        T = min(T, this->tri11(VAL(l0), VAL(l1), s_hat, h));
      }
      if (nb[l2]) {
        s_hat = this->estimate_speed(s, s_[l0], s_[l2]);
        T = min(T, this->tri11(VAL(l0), VAL(l2), s_hat, h));
      }
      if (nb[l1] && nb[l2]) {
        s_hat = this->estimate_speed(s, s_[l1], s_[l2]);
        T = min(T, this->tetra111(VAL(l0), VAL(l1), VAL(l2), s_hat, h));
      }
    }
  }
  if (nb[0] && nb[2] && nb[4]) {
    s_hat = this->estimate_speed(s, s_[0], s_[2], s_[4]);
    T = min(T, this->tetra111(VAL(0), VAL(2), VAL(4), s_hat, h));
  }
  if (nb[1] && nb[3] && nb[5]) {
    s_hat = this->estimate_speed(s, s_[1], s_[3], s_[5]);
    T = min(T, this->tetra111(VAL(1), VAL(3), VAL(5), s_hat, h));
  }

#ifdef PRINT_UPDATES
  printf("olim6_rect::update_impl: T <- %g\n", T);
#endif
}

#endif // __OLIM6_RECT_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
