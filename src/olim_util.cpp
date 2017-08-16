#include "olim_util.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

static void check_params(double u0, double u1, double h) {
  (void) u0;
  (void) u1;
  (void) h;
  assert(h > 0);
  assert(u0 >= 0);
  assert(u1 >= 0);
  assert(!std::isinf(u0));
  assert(!std::isinf(u1));
}

/**
 * Adjacent triangle update with constant quadrature rule (used by
 * olim8_rhr and olim8_mp0c).
 */
double rhr_adj(double u0, double u1, double s_est, double h, double * lam) {
  check_params(u0, u1, h);
  assert(s_est >= 0);

  double c = (u0 - u1)/(s_est*h);
  // TODO: use copysign for next line instead?
  double _lam = 0.5 + (c > 0 ? 1 : -1)*std::fabs(c)/(2*std::sqrt(2 - c*c));
  if (lam != nullptr) {
    *lam = _lam;
  }
  return (1 - _lam)*u0 + _lam*u1 + s_est*h*sqrt(2*_lam*(_lam - 1) + 1);
}

/**
 * Diagonal triangle update with constant quadrature rule (used by
 * olim8_rhr and olim8_mp0c).
 */
double rhr_diag(double u0, double u1, double s_est, double h) {
  check_params(u0, u1, h);
  assert(s_est >= 0);
  
  // TODO: make this one as simple as the adj rule
  double c = std::fabs(u0 - u1)/(s_est*h);
  double sgn = u0 > u1 ? 1 : -1;
  double lam = std::max(0.0, std::min(1.0, sgn*c/std::sqrt(1 - c*c)));
  return (1 - lam)*u0 + lam*u1 + s_est*h*sqrt(lam*lam + 1);
}

static double polyval(double * coefs, int ncoefs, double x) {
  double y = 0;
  for (int i = ncoefs - 1; i > 0; --i) {
    y += coefs[i];
    y *= x;
  }
  return y + coefs[0];
}

static int sigma(double ** polys, double x) {
  // TODO: we can optimize this by packing bits probably
  // int nsigns = 0;
  // int signs[5];
  // double y;
  // for (int i = 0; i < 5; ++i) {
  //   y = polyval(polys[i], 5 - i, x);
  //   if (y > 0) {
  //     signs[nsigns++] = 1;
  //   } else if (y < 0) {
  //     signs[nsigns++] = -1;
  //   }
  // }
  // int changes = 0;
  // for (int i = 1, j = 0; i < nsigns; ++i, ++j) {
  //   if (signs[i] != signs[j]) {
  //     ++changes;
  //   }
  // }
  // return changes;

  char signs = polyval(polys[0], 5, x) > 0 ? 127 : 0;

  double y;
  for (int i = 1; i < 5; ++i) {
    y = polyval(polys[i], 5 - i, x);
    if (y > 0) {
      signs = (signs << 1) | 1;
    } else if (y < 0) {
      signs <<= 1;
    }
  }

  return ((signs & 1) ^ (signs & 2)) +
    ((signs & 2) ^ (signs & 3)) +
    ((signs & 3) ^ (signs & 4)) +
    ((signs & 4) ^ (signs & 5));
}

int sturm(double ** polys, double l, double r) {
  int nroots = 0, degree = 4;
  double * a = polys[0];
  while (a[degree] < 1e-13) {
    --degree;
  }
  assert(degree > 1);
  if (degree == 2) {
    double disc = a[1]*a[1] - 4*a[2]*a[0];
    if (disc < 0) {
      nroots = 0;
    } else if (disc == 0) {
      double x0 = (-a[1] + std::sqrt(disc))/(2*a[2]);
      nroots = l < x0 && x0 <= r ? 1 : 0;
    } else {
      double x0 = (-a[1] + std::sqrt(disc))/(2*a[2]),
        x1 = (-a[1] - std::sqrt(disc))/(2*a[2]);
      nroots = (int) (l < x0 && x0 <= r) + (int) (l < x1 && x1 <= r);
    }
  }
  assert(degree != 3);
  if (degree == 4) {
    nroots = sigma(polys, l) - sigma(polys, r);
  }
  return nroots;
}

static double secant(double * poly, double x0, double x1, double l, double r,
                     bool & foundroot, double tol = 1e-13) {
  double x, f0, f1;
  do {
    f0 = polyval(poly, 5, x0);
    f1 = polyval(poly, 5, x1);
    x = (x1*f0 - x0*f1)/(f0 - f1);
    if (x < l || r < x) {
      foundroot = false;
      return x;
    }
    x1 = x0;
    x0 = x;
  } while (fabs(polyval(poly, 5, x)) > tol);
  foundroot = true;
  return x;
}

static void rec(double ** polys, double * roots, int & root,
                double l, double r) {
  int nroots = sturm(polys, l, r);
  if (nroots == 1) {
    double h = 0.1, x0;
    bool foundroot = false;
    x0 = secant(polys[0], l, l + (r - l)*h, l, r, foundroot);
    if (!foundroot) {
      x0 = secant(polys[0], r, r + (l - r)*h, l, r, foundroot);
    }
    assert(foundroot);
    roots[root++] = x0;
  } else if (nroots > 1) {
    double mid = (l + r)/2;
    rec(polys, roots, root, l, mid);
    rec(polys, roots, root, mid, r);
  }
}

void find_quartic_roots(double * a, double * roots, double l, double r) {
  // TODO: use descartes sign rule here to save flops?
  
  double b[4] = {a[1], 2*a[2], 3*a[3], 4*a[4]};
  double c[3] = {
    -a[0] + a[1]*a[3]/(16*a[4]),
    (-6*a[1] + a[2]*a[3]/a[4])/8,
    -a[2]/2 + 3*a[3]*a[3]/(16*a[4])
  };
  double d[2] = {
    -b[0] + c[0]*(-b[3]*c[1] + b[2]*c[2])/(c[2]*c[2]),
    (-b[3]*(c[1]*c[1] + c[0]*c[2]) + (b[2]*c[1] - b[1]*c[2])*c[2])/(c[2]*c[2])
  };
  double e[1] = {-c[0] + d[0]*(-c[2]*d[0] + c[1]*d[1])/(d[1]*d[1])};

  // printf("double a[5] = {%g, %g, %g, %g, %g}\n", a[0], a[1], a[2], a[3], a[4]);
  // printf("double b[4] = {%g, %g, %g, %g}\n", b[0], b[1], b[2], b[3]);
  // printf("double c[3] = {%g, %g, %g}\n", c[0], c[1], c[2]);
  // printf("double d[2] = {%g, %g}\n", d[0], d[1]);
  // printf("double e[1] = {%g}\n", e[0]);
  // printf("\n");

  double * polys[5] = {a, b, c, d, e};

  // TODO: check if f(0) = 0? (since sturm handles the half-open
  // interval (a, b])

  int root = 0;
  rec(polys, roots, root, l, r);

  // printf("roots = {");
  // if (root > 0) {
  //   for (int i = 0; i < root - 1; ++i) {
  //     printf("%g, ", roots[i]);
  //   }
  //   printf("%g}\n", roots[root - 1]);
  // } else {
  //   printf("}\n");
  // }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
