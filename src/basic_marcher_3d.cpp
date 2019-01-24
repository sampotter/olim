#include "basic_marcher_3d.hpp"

#include <assert.h>
#include <math.h>

#include "common.hpp"

// TODO: remove these macros and replace with lambdas

#define COMPUTE_DISC_2PT() (2*sh_sq - (T1 - T2)*(T1 - T2))

#define COMPUTE_VALUE_2PT() ((T1 + T2 + sqrt(disc))/2)

#define COMPUTE_DISC_3PT() \
  (3*sh_sq - 2*(T1*T1 + T2*T2 + T3*T3 - T1*T2 - T1*T3 - T2*T3))

#define COMPUTE_VALUE_3PT() ((T1 + T2 + T3 + sqrt(disc))/3)

void basic_marcher_3d::update_impl(int lin_hat, int * nb, int parent, double & T)
{
  (void) lin_hat;
  (void) parent;

  // int i = n->get_i(), j = n->get_j(), k = n->get_k();

  double sh = get_h()*s_hat, sh_sq = sh*sh;
  double T1 = 0, T2 = 0, T3 = 0, disc = 0;

  auto const value = [&] (int i) {
    return _U[nb[i]];
  };

  for (int l0 = 0, l1 = 1, l2 = 2; l0 < 6;
       ++l0, l1 = (l1 + 1) % 6, l2 = (l2 + 1) % 6) {
    if (nb[l0] != -1) {
      T1 = value(l0);
      T = fmin(T, T1 + sh);
      if (nb[l1] != -1) {
        T2 = value(l1);
        disc = COMPUTE_DISC_2PT();
        if (disc > 0) T = fmin(T, COMPUTE_VALUE_2PT());
      }
      if (nb[l2] != -1) {
        T2 = value(l2);
        disc = COMPUTE_DISC_2PT();
        if (disc > 0) T = fmin(T, COMPUTE_VALUE_2PT());
      }
      if (nb[l1] != -1 && nb[l2] != -1) {
        T2 = value(l1), T3 = value(l2);
        disc = COMPUTE_DISC_3PT();
        if (disc > 0) T = fmin(T, COMPUTE_VALUE_3PT());
      }
    }
  }
  if (nb[0] != -1 && nb[2] != -1 && nb[4] != -1) {
    T1 = value(0), T2 = value(2), T3 = value(4);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) T = fmin(T, COMPUTE_VALUE_3PT());
  }
  if (nb[1] != -1 && nb[3] != -1 && nb[5] != -1) {
    T1 = value(1), T2 = value(3), T3 = value(5);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) T = fmin(T, COMPUTE_VALUE_3PT());
  }
}
