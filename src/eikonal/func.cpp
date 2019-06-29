#include "func.hpp"

#include <assert.h>

// TODO: we only have to compute Lagrange multipliers at interior
// points... we can try to optimize this accordingly

template <>
void lagmults<2>(vec2<double> const & lam, vec2<double> const & df,
                 double const * d2f, vec2<double> & mu, int * k_ptr)
{
  check_lambda<2>(lam);

  int I[2] = {-1, -1};

  int & k = *k_ptr;

  // First, find the active set, I, for lam
  k = 0;
  if (fabs(lam[0]) < eps<double>) I[k++] = 0;
  if (fabs(lam[1]) < eps<double>) I[k++] = 1;
  if (fabs(1 - lam[0] - lam[1]) < eps<double>) I[k++] = 2;

  // TODO: we could try just adding an assert here that k < 2

  if (k == 2 && I[0] != 2 && I[1] != 2) {
    // If I = {0, 1}, then mu = df. Simple enough.
    mu[0] = df[0];
    mu[1] = df[1];
  } else {
    if (k == 1) {
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
    } else if (k == 2) {
      if (fabs(d2f[0]*d2f[2] - d2f[1]*d2f[1]) < eps<double>) {
        mu[0] = mu[1] = 0;
      } else {
        if (I[0] != 1 && I[1] != 1) {
          // When I = {0, 2}, A = [1 0; -1 -1] = inv(A), which we use to
          // compute the following matrix-vector product.
          mu[0] = df[0] - df[1];
          mu[1] = -df[1];
        } else if (I[0] != 0 && I[1] != 0) {
          // Likewise, when I = {1, 2}, we have A = [0 1; -1 -1] and
          // inv(A) = [-1 -1; 1 0], which we use to compute mu.
          mu[0] = df[1] - df[0];
          mu[1] = -df[0];
        } else {
          assert(false);
        }
      }
    } else {
      assert(false);
    }
  }

#if OLIM_DEBUG && !RELWITHDEBINFO
  assert(!isnan(mu[0]));
  assert(!isnan(mu[1]));
#endif
}
