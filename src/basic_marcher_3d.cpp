#include "basic_marcher_3d.hpp"

#include <assert.h>
#include <math.h>

#include "common.hpp"

// TODO: remove these macros and replace with lambdas

#define COMPUTE_DISC_2PT() (2*sh_sq - (U1 - U2)*(U1 - U2))

#define COMPUTE_VALUE_2PT() ((U1 + U2 + sqrt(disc))/2)

#define COMPUTE_DISC_3PT() \
  (3*sh_sq - 2*(U1*U1 + U2*U2 + U3*U3 - U1*U2 - U1*U3 - U2*U3))

#define COMPUTE_VALUE_3PT() ((U1 + U2 + U3 + sqrt(disc))/3)

void basic_marcher_3d::update_impl(int lin_hat, int * nb, int parent, double & U)
{
  (void) lin_hat;
  (void) parent;

  double sh = get_h()*s_hat, sh_sq = sh*sh;
  double U1 = 0, U2 = 0, U3 = 0, disc = 0;

  auto const value = [&] (int i) {
    return _U[nb[i]];
  };

  for (int l0 = 0, l1 = 1, l2 = 2; l0 < 6;
       ++l0, l1 = (l1 + 1) % 6, l2 = (l2 + 1) % 6) {
    if (nb[l0] != -1) {
      U1 = value(l0);
      U = fmin(U, U1 + sh);
      if (nb[l1] != -1) {
        U2 = value(l1);
        disc = COMPUTE_DISC_2PT();
        if (disc > 0) U = fmin(U, COMPUTE_VALUE_2PT());
      }
      if (nb[l2] != -1) {
        U2 = value(l2);
        disc = COMPUTE_DISC_2PT();
        if (disc > 0) U = fmin(U, COMPUTE_VALUE_2PT());
      }
      if (nb[l1] != -1 && nb[l2] != -1) {
        U2 = value(l1), U3 = value(l2);
        disc = COMPUTE_DISC_3PT();
        if (disc > 0) U = fmin(U, COMPUTE_VALUE_3PT());
      }
    }
  }
  if (nb[0] != -1 && nb[2] != -1 && nb[4] != -1) {
    U1 = value(0), U2 = value(2), U3 = value(4);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) U = fmin(U, COMPUTE_VALUE_3PT());
  }
  if (nb[1] != -1 && nb[3] != -1 && nb[5] != -1) {
    U1 = value(1), U2 = value(3), U3 = value(5);
    disc = COMPUTE_DISC_3PT();
    if (disc > 0) U = fmin(U, COMPUTE_VALUE_3PT());
  }
}
