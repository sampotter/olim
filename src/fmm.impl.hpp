#pragma once

#include <assert.h>
#include <math.h>

#include "common.hpp"

template <>
void fmm<2>::update_impl(int lin, int const * nb, int parent, double & U) {
  (void) lin;
  (void) parent;

  double sh = this->_s[lin]*this->get_h();

  double U1 = fmin(
    nb[0] != -1 ? _U[nb[0]] : inf<double>,
    nb[2] != -1 ? _U[nb[2]] : inf<double>);

  double U2 = fmin(
    nb[1] != -1 ? _U[nb[1]] : inf<double>,
    nb[3] != -1 ? _U[nb[3]] : inf<double>);

  bool U1_inf = isinf(U1), U2_inf = isinf(U2);

  if (!U1_inf && !U2_inf) {
    double diff = U1 - U2, disc = 2*sh*sh - diff*diff;
    U = disc > 0 ? fmin(U, (U1 + U2 + sqrt(disc))/2) : U;
  } else if (U1_inf) {
    U = fmin(U, U2 + sh);
  } else if (U2_inf) {
    U = fmin(U, U1 + sh);
  } else {
    assert(false);
  }
}

template <>
void fmm<3>::update_impl(int lin_hat, int const * nb, int parent, double & U)
{
  (void) lin_hat;
  (void) parent;

  double sh = get_h()*this->_s[lin_hat], sh_sq = sh*sh;
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

  int l1s[7] = {1, 2, 3, 4, 5, 0, 1};
  int * l2s = &l1s[1];
  for (int l0 = 0, l1, l2; l0 < 6; ++l0) {
    if (nb[l0] != -1) {
      l1 = l1s[l0];
      l2 = l2s[l0];
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
