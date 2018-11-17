#include "cost_funcs.hpp"

#include <cassert>
#include <cmath>
#include <utility>

#define __invert2x2inplace(X) do {              \
    double const det = X[0]*X[2] - X[1]*X[1];   \
    std::swap(X[0], X[2]);                      \
    X[0] /= det;                                \
    X[1] /= -det;                               \
    X[2] /= det;                                \
  } while (0)

template <>
void lagmults<2>(double const * lam, double const * df, double const * d2f,
                 double * mu, int * k)
{
  // First, find the active set for lam

  int active_set[2] = {-1, -1};
  int num_active = 0;

  if (fabs(lam[1]) < EPS(double)) active_set[num_active++] = 1;
  if (fabs(lam[0]) < EPS(double)) active_set[num_active++] = 0;
  if (fabs(1 - lam[0] - lam[1]) < EPS(double)) active_set[num_active++] = 2;

  *k = num_active;

  // If the active_set is {0, 1}, then mu = -df, which means that we
  // can return early.

  if (num_active == 2 &&
      ((active_set[0] == 0 && active_set[1] == 1) ||
       (active_set[0] == 1 && active_set[1] == 0))) {
    mu[0] = -df[0];
    mu[1] = -df[1];
  } else {
    double inv_d2f[3] = {d2f[0], d2f[1], d2f[2]};
    __invert2x2inplace(inv_d2f);

    double p[2];
    p[0] = inv_d2f[0]*df[0] + inv_d2f[1]*df[1];
    p[1] = inv_d2f[1]*df[0] + inv_d2f[2]*df[1];

    if (num_active == 1) {
      if (active_set[0] == 0) {
        mu[0] = -p[0]/d2f[0];
      } else if (active_set[0] == 1) {
        mu[0] = -p[1]/d2f[2];
      } else if (active_set[0] == 2) {
        mu[0] = (p[0] + p[1])/(d2f[0] + 2*d2f[1] + d2f[2]);
      } else {
        assert(false);
        mu[0] = mu[1] = INF(double);
      }
    } else if (num_active == 2) {
      // TODO: we can do this in place without defining A and b
      double A[3], b[2];
      if ((active_set[0] == 0 && active_set[1] == 2) ||
          (active_set[0] == 2 && active_set[1] == 0)) {
        A[0] = d2f[0];
        A[1] = -(d2f[0] + d2f[1]);
        A[2] = d2f[0] + 2*d2f[1] + d2f[2];
        __invert2x2inplace(A);
        b[0] = -p[0];
        b[1] = p[0] + p[1];
      } else if ((active_set[0] == 1 && active_set[1] == 2) ||
                 (active_set[0] == 2 && active_set[1] == 1)) {
        A[0] = d2f[2];
        A[1] = -(d2f[1] + d2f[2]);
        A[2] = d2f[0] + 2*d2f[1] + d2f[2];
        __invert2x2inplace(A);
        b[0] = -p[1];
        b[1] = p[0] + p[1];
      } else {
        assert(false);
        A[0] = A[1] = A[2] = b[0] = b[1] = INF(double);
      }
      mu[0] = A[0]*b[0] + A[1]*b[1];
      mu[1] = A[1]*b[0] + A[2]*b[1];
    } else {
      assert(false);
      mu[0] = mu[1] = INF(double);
    }
  }

  CHECK(mu[0]);
  CHECK(mu[1]);
}

#undef __invert2x2inplace
