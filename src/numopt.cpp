#include "numopt.hpp"

#include <assert.h>
#include <math.h>

template <>
void
qpe_bary<2, 0>(double const * G, double const * c, double * x)
{
  x[0] = 0.0;
  x[1] = -c[1]/G[2];
}

template <>
void
qpe_bary<2, 1>(double const * G, double const * c, double * x)
{
  x[0] = -c[0]/G[0];
  x[1] = 0.0;
}

template <>
void
qpe_bary<2, 2>(double const * G, double const * c, double * x)
{
  double pZ = ((G[2] - G[0])/2 + c[1] - c[0])/(G[0] - 2*G[1] + G[2]);
  x[0] = 0.5 + pZ;
  x[1] = 0.5 - pZ;
}

#define __compute_p() do {                          \
    double det = G[0]*G[2] - G[1]*G[1];             \
    p[0] = -(x[0] + (G[2]*c[0] - G[1]*c[1])/det);   \
    p[1] = -(x[1] + (G[0]*c[1] - G[1]*c[0])/det);   \
  } while (0)

#define __compute_y() do {                      \
    y[0] = G[0]*x[0] + G[1]*x[1] + c[0];        \
    y[1] = G[1]*x[0] + G[2]*x[1] + c[1];        \
  } while (0)

#define __num_active()                                              \
  ((active[0] ? 1 : 0) + (active[1] ? 1 : 0) + (active[2] ? 1 : 0))

template <>
void
qpi_bary<2>(double const * G, double const * c, double const * x0,
            double * x, bool * error, double tol, int niters)
{
  assert(x != nullptr);

  if (error) {
    *error = false;
  }

  double xprev[2], p[2] = {0, 0}, y[2], alpha, alpha_new;
  if (x0) {
    x[0] = x0[0];
    x[1] = x0[1];
  } else {
    x[0] = x[1] = 0.0;
  }

  bool active[3], ind[3];
  active[0] = fabs(x[0]) <= tol;
  active[1] = fabs(x[1]) <= tol;
  active[2] = fabs(1 - x[0] - x[1]) <= tol;

  int k = 0, num_active, argmin;
  while (true) {
    xprev[0] = x[0];
    xprev[1] = x[1];
    num_active = __num_active();

    __compute_y();

    if (num_active == 0) {
      __compute_p();
    } else if (num_active == 1) {
      if (active[0]) qpe_bary<2, 0>(G, c, p);
      else if (active[1]) qpe_bary<2, 1>(G, c, p);
      else if (active[2]) qpe_bary<2, 2>(G, c, p);
      else assert(false);
      p[0] -= xprev[0];
      p[1] -= xprev[1];
    } else if (num_active == 2) {
      p[0] = p[1] = 0.0;
    }

    if (fmax(fabs(p[0]), fabs(p[1])) <= tol) {
      if (num_active == 0) {
        break;
      } else if (num_active == 1) {
        if (active[0]) {
          if (y[0] >= 0) break;
          else active[0] = false;
        } else if (active[1]) {
          if (y[1] >= 0) break;
          else active[1] = false;
        } else if (active[2]) {
          if (y[0] + y[1] <= 0) break;
          else active[2] = false;
        } else {
          assert(false);
        }
      } else if (num_active == 2) {
        if (active[0] && active[1]) {
          if (y[0] >= 0 && y[1] >= 0) break;
          else active[y[0] < y[1] ? 0 : 1] = false;
        } else if (active[0] && active[2]) {
          if (y[1] <= 0 && y[1] <= y[0]) break;
          else active[y[0] < 0 ? 0 : 2] = false;
        } else if (active[1] && active[2]) {
          if (y[0] <= 0 && y[0] <= y[1]) break;
          else active[y[1] < 0 ? 1 : 2] = false;
        } else {
          assert(false);
        }
      } else {
        assert(false);
      }
    } else {
      alpha = 1;
      argmin = -1;
      ind[0] = !active[0] && p[0] < 0;
      ind[1] = !active[1] && p[1] < 0;
      ind[2] = !active[2] && p[0] + p[1] > 0;
      if (ind[0] || ind[1] || ind[2]) {
        if (ind[0]) {
          alpha_new = -x[0]/p[0];
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 0;
          }
        }
        if (ind[1]) {
          alpha_new = -x[1]/p[1];
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 1;
          }
        }
        if (ind[2]) {
          alpha_new = (1 - x[0] - x[1])/(p[0] + p[1]);
          if (alpha_new < alpha) {
            alpha = alpha_new;
            argmin = 2;
          }
        }
        alpha = fmax(0.0, fmin(alpha, 1.0));
        if (alpha < 1) {
          active[argmin] = true;
        }
      }
      x[0] += alpha*p[0];
      x[1] += alpha*p[1];
    }

    if (++k == niters) {
      if (error) *error = true;
      break;
    }
  }
}

#undef __compute_p
#undef __compute_y
#undef __num_active
