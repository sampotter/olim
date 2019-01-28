#include "basic_marcher_3d.hpp"

#include <assert.h>
#include <math.h>

#include "common.hpp"

void basic_marcher_3d::update_impl(int lin_hat, int * nb, int parent, double & U)
{
  (void) lin_hat;
  (void) parent;

  double sh = get_h()*s_hat, sh_sq = sh*sh;
  double U1 = 0, U2 = 0, U3 = 0, disc = 0;

  auto const value = [&] (int i) {
    return _U[nb[i]];
  };

  auto const disc_2pt = [&] () {
    return 2*sh_sq - (U1 - U2)*(U1 - U2);
  };

  auto const disc_3pt = [&] () {
    return 3*sh_sq - 2*(U1*U1 + U2*U2 + U3*U3 - U1*U2 - U1*U3 - U2*U3);
  };

  auto const value_2pt = [&] () {
    return (U1 + U2 + sqrt(disc))/2;
  };

  auto const value_3pt = [&] () {
    return (U1 + U2 + U3 + sqrt(disc))/3;
  };

  for (int l0 = 0, l1 = 1, l2 = 2; l0 < 6;
       ++l0, l1 = (l1 + 1) % 6, l2 = (l2 + 1) % 6) {
    if (nb[l0] != -1) {
      U1 = value(l0);
      U = fmin(U, U1 + sh);
      if (nb[l1] != -1) {
        U2 = value(l1);
        disc = disc_2pt();
        if (disc > 0) {
          U = fmin(U, value_2pt());
        }
      }
      if (nb[l2] != -1) {
        U2 = value(l2);
        disc = disc_2pt();
        if (disc > 0) {
          U = fmin(U, value_2pt());
        }
      }
      if (nb[l1] != -1 && nb[l2] != -1) {
        U2 = value(l1), U3 = value(l2);
        disc = disc_3pt();
        if (disc > 0) {
          U = fmin(U, value_3pt());
        }
      }
    }
  }
  if (nb[0] != -1 && nb[2] != -1 && nb[4] != -1) {
    U1 = value(0), U2 = value(2), U3 = value(4);
    disc = disc_3pt();
    if (disc > 0) {
      U = fmin(U, value_3pt());
    }
  }
  if (nb[1] != -1 && nb[3] != -1 && nb[5] != -1) {
    U1 = value(1), U2 = value(3), U3 = value(5);
    disc = disc_3pt();
    if (disc > 0) {
      U = fmin(U, value_3pt());
    }
  }
}
