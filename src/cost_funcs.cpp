#include "cost_funcs.hpp"

#include <assert.h>

template <>
void lagmults<2>(double const * lam, double const * df, double const * d2f,
                 double * mu, int * k)
{
  int I[2] = {-1, -1};

  // First, find the active set, I, for lam
  *k = 0;
  if (fabs(lam[0]) < eps<double>) I[(*k)++] = 0;
  if (fabs(lam[1]) < eps<double>) I[(*k)++] = 1;
  if (fabs(1 - lam[0] - lam[1]) < eps<double>) I[(*k)++] = 2;

  if (*k == 2 && I[0] != 2 && I[1] != 2) {
    // If I = {0, 1}, then mu = df. Simple enough.
    mu[0] = df[0];
    mu[1] = df[1];
  } else {
    if (*k == 1) {
      // When k = 1, things simplify a lot:
      if (I[0] == 0) {
        mu[0] = df[0] - d2f[1]*df[1]/d2f[2];
      } else if (I[0] == 1) {
        mu[0] = df[1] - d2f[1]*df[0]/d2f[0];
      } else if (I[0] == 2) {
        double A = d2f[2] - d2f[1], B = d2f[0] - d2f[1];
        mu[0] = -(A*df[0] + B*df[1])/(A + B);
      } else {
        assert(false);
      }
    } else if (*k == 2) {
      double det, inv[3], p[2];
      det = d2f[0]*d2f[2] - d2f[1]*d2f[1];
      if (fabs(det) < eps<double>) {
        mu[0] = mu[1] = 0;
      } else {
        inv[0] = d2f[2]/det;
        inv[1] = -d2f[1]/det;
        inv[2] = d2f[0]/det;
        p[0] = inv[0]*df[0] + inv[1]*df[1];
        p[1] = inv[1]*df[0] + inv[2]*df[1];
        if (I[0] != 1 && I[1] != 1) {
          // When I = {0, 2}, A = [1 0; -1 -1] = inv(A), which we use to
          // compute the following matrix-vector product.
          mu[0] = (d2f[0] - 2*d2f[1] + d2f[2])*p[0] + (d2f[1] - d2f[2])*(p[0] + p[1]);
          mu[1] = (d2f[2] - d2f[1])*p[0] - d2f[2]*(p[0] + p[1]);
        } else if (I[0] != 0 && I[1] != 0) {
          // Likewise, when I = {1, 2}, we have A = [0 1; -1 -1] and
          // inv(A) = [-1 -1; 1 0], which we use to compute mu.
          mu[0] = (d2f[0] - 2*d2f[1] + d2f[2])*p[1] - (d2f[0] - d2f[1])*(p[0] + p[1]);
          mu[1] = (d2f[0] - d2f[1])*p[1] - d2f[0]*(p[0] + p[1]);
        } else {
          assert(false);
        }
      }
    } else {
      assert(false);
    }
  }

  assert(!isnan(mu[0]));
  assert(!isnan(mu[1]));
}
