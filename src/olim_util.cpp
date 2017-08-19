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

static char nbits[16] = {
  0, // 0000
  1, // 0001
  1, // 0010
  2, // 0011
  1, // 0100
  2, // 0101
  2, // 0110
  3, // 0111
  1, // 1000
  2, // 1001
  2, // 1010
  3, // 1011
  2, // 1100
  3, // 1101
  3, // 1110
  4  // 1111
};

int sigma(double ** polys, double x) {
  double * a = polys[0];
  double y = a[0] + x*(a[1] + x*(a[2] + x*(a[3] + x*a[4])));
  char signs = y > 0 ? 127 : 0;

  a = polys[1];
  y = a[0] + x*(a[1] + x*(a[2] + x*a[3]));
  signs = (signs << (y != 0)) | (y > 0);

  a = polys[2];
  y = a[0] + x*(a[1] + x*a[2]);
  signs = (signs << (y != 0)) | (y > 0);

  a = polys[3];
  y = a[0] + x*a[1];
  signs = (signs << (y != 0)) | (y > 0);

  a = polys[4];
  y = a[0];
  signs = (signs << (y != 0)) | (y > 0);

  // TODO: compare with modding and adding---that may actually be
  // faster
  return nbits[((signs ^ (signs << 1)) >> 1) & 15];
}

int oldsigma(double ** polys, double x) {
  int nsigns = 0;
  int signs[5];
  double y;
  for (int i = 0; i < 5; ++i) {
    y = polyval(polys[i], 5 - i, x);
    if (y > 0) {
      signs[nsigns++] = 1;
    } else if (y < 0) {
      signs[nsigns++] = -1;
    }
  }
  int changes = 0;
  for (int i = 1, j = 0; i < nsigns; ++i, ++j) {
    if (signs[i] != signs[j]) {
      ++changes;
    }
  }
  return changes;
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
      // TODO: use unsigned form for comparison
      nroots = l < x0 && x0 <= r ? 1 : 0;
    } else {
      double x0 = (-a[1] + std::sqrt(disc))/(2*a[2]),
        x1 = (-a[1] - std::sqrt(disc))/(2*a[2]);
      // TODO: use unsigned form for comparison
      nroots = (int) (l < x0 && x0 <= r) + (int) (l < x1 && x1 <= r);
    }
  }
  assert(degree != 3);
  if (degree == 4) {
    // TODO: Sturm's theorem finds roots on (l, r]; the case p(l) = 0
    // *does* happen and *isn't* handled by our fast sigma---so,
    // perturb l to the right a very small amount to get around this
    double l_plus_eps = l + std::numeric_limits<decltype(l)>::epsilon();
    assert(sigma(polys, l_plus_eps) == oldsigma(polys, l_plus_eps));
    assert(sigma(polys, r) == oldsigma(polys, r));
    nroots = sigma(polys, l_plus_eps) - sigma(polys, r);
  }
  return nroots;
}

// TODO: remove use of polyval

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

static void findroot(double * a, double l, double r, double h, double * x0) {
  bool found = false;
  *x0 = secant(a, l, l + (r - l)*h, l, r, found);
  if (!found) *x0 = secant(a, r, r + (l - r)*h, l, r, found);
  assert(found);
}

struct interval {
  interval() {}
  interval(int nroots, double l, double r): nroots {nroots}, l {l}, r {r} {}
  int nroots {-1};
  double l {-1}, r {-1};
};

void find_quartic_roots(double * a, double * roots, double l, double r) {
  // TODO: simplify arithmetic (low priority)

  double b[4] = {a[1], 2*a[2], 3*a[3], 4*a[4]}; // 1st der. of a
  double c[3] = { // -rem(a, b)
    -a[0] + a[1]*a[3]/(16*a[4]),
    (-6*a[1] + a[2]*a[3]/a[4])/8,
    -a[2]/2 + 3*a[3]*a[3]/(16*a[4])
  };
  double d[2] = { // -rem(b, c)
    -b[0] + c[0]*(-b[3]*c[1] + b[2]*c[2])/(c[2]*c[2]),
    (-b[3]*(c[1]*c[1] + c[0]*c[2]) + (b[2]*c[1] - b[1]*c[2])*c[2])/(c[2]*c[2])
  };
  double e[1] = { // -rem(c, d)
    -c[0] + d[0]*(-c[2]*d[0] + c[1]*d[1])/(d[1]*d[1])
  };
  double f[3] = {b[1], 2*b[2], 3*b[3]}; // 2nd der. of a

  double * polys[6] = {a, b, c, d, e, f};

  int root = 0;
  
  if (fabs(a[0]) < 1e-14) {
    roots[root++] = 0;
  }

  double mid;
  int nroots = sturm(polys, l, r);
  interval ivals[4];
  int i = 0;
  if (nroots >= 1) {
    ivals[i++] = {nroots, l, r};
  }
  while (i > 0) {
    interval ival = ivals[--i];
    if (ival.nroots == 1) {
      findroot(polys[0], ival.l, ival.r, 0.1, &roots[root++]);
    } else if (ival.nroots > 1) {
      mid = (l + r)/2;
      int sigl = sigma(polys, l), sigm = sigma(polys, mid),
        sigr = sigma(polys, r);
      if ((nroots = sturm(polys, l, mid)) >= 1) {
        assert(0 <= i && i < 4);
        ivals[i++] = {nroots, l, mid};
      }
      if ((nroots = sturm(polys, mid, r)) >= 1) {
        assert(0 <= i && i < 4);
        ivals[i++] = {nroots, mid, r};
      }
    }
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
