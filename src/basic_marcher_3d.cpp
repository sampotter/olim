#include "basic_marcher_3d.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "common.macros.hpp"

#define COMPUTE_DISC_2PT() (2*sh*sh - (T1 - T2)*(T1 - T2))

#define COMPUTE_VALUE_2PT() ((T1 + T2 + std::sqrt(disc))/2)

#define COMPUTE_DISC_3PT() \
  (3*sh_sq - 2*(T1*T1 + T2*T2 + T3*T3 - T1*T2 - T1*T3 - T2*T3))

#define COMPUTE_VALUE_3PT() ((T1 + T2 + T3 + std::sqrt(disc))/3)

void basic_marcher_3d::update_impl(
  node_3d * n, node_3d ** nb, int parent, double & T)
{
  (void) parent;

  using std::min;

  int i = n->get_i(), j = n->get_j(), k = n->get_k();

  double sh = get_h()*get_speed(i, j, k), sh_sq = sh*sh;
  double T1 = 0, T2 = 0, T3 = 0, disc = 0;

  for (int l0 = 0, l1 = 1, l2 = 2; l0 < 6;
       ++l0, l1 = (l1 + 1) % 6, l2 = (l2 + 1) % 6) {
    if (nb[l0]) {
      T1 = VAL(l0);
      T = min(T, T1 + sh);
      if (nb[l1]) {
        T2 = VAL(l1);
        disc = COMPUTE_DISC_2PT();
        if (disc > 0) T = min(T, COMPUTE_VALUE_2PT());
      }
      if (nb[l2]) {
        T2 = VAL(l2);
        disc = COMPUTE_DISC_2PT();
        if (disc > 0) T = min(T, COMPUTE_VALUE_2PT());
      }
      if (nb[l1] && nb[l2]) {
        T2 = VAL(l1), T3 = VAL(l2);
        disc = COMPUTE_DISC_3PT();
        if (disc > 0) T = min(T, COMPUTE_VALUE_3PT());
      }
    }
  }
  if (nb[0] && nb[2] && nb[4]) {
    T1 = VAL(0), T2 = VAL(2), T3 = VAL(4);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) T = min(T, COMPUTE_VALUE_3PT());
  }
  if (nb[1] && nb[3] && nb[5]) {
    T1 = VAL(1), T2 = VAL(3), T3 = VAL(5);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) T = min(T, COMPUTE_VALUE_3PT());
  }
}
