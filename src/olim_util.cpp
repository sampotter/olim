#include "olim_util.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

#include <gsl/gsl_poly.h>

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

double mp0l_adj(double u0, double u1, double s, double s0, double s1, double h) {
  check_params(u0, u1, h);
  assert(s >= 0);
  assert(s0 >= 0);
  assert(s1 >= 0);

  double s0bar = (s + s0)/2;
  double ds = s1 - s0;

  if (s == s0 && s == s1) {
    return rhr_adj(u0, u1, s, h);
  } else if (ds == 0) {
    return rhr_adj(u0, u1, s0bar, h);
  }

  double s1bar = (s + s1)/2;
  double du = u1 - u0;
  double b = s0bar/ds - 0.75;
  double c = (ds - s - s0)/(4*ds) + du/(2*ds*h);
  double disc = b*b - 4*c;
  assert(disc >= 0);
  double lam1 = (-b + std::sqrt(disc))/2;
  double lam2 = (-b - std::sqrt(disc))/2;

  double T = std::numeric_limits<double>::infinity(), T_new, one_minus_lam,
    lam_sq, dist, slambar;
  if (0 <= lam1 && lam1 <= 1) {
    one_minus_lam = 1 - lam1;
    lam_sq = lam1*lam1;
    slambar = (one_minus_lam*s0bar + lam1*s1bar);
    dist = h*std::sqrt(lam_sq + one_minus_lam*one_minus_lam);
    T_new = one_minus_lam*u0 + lam1*u1 + slambar*dist;
    T = std::min(T, T_new);
  }
  if (0 <= lam2 && lam2 <= 1) {
    one_minus_lam = 1 - lam2;
    lam_sq = lam2*lam2;
    slambar = (one_minus_lam*s0bar + lam2*s1bar);
    dist = h*std::sqrt(lam_sq + one_minus_lam*one_minus_lam);
    T_new = one_minus_lam*u0 + lam2*u1 + slambar*dist;
    T = std::min(T, T_new);
  }
  return T;
}

double mp0l_diag(double u0, double u1, double s, double s0, double s1, double h) {
  check_params(u0, u1, h);
  assert(s >= 0);
  assert(s0 >= 0);
  assert(s1 >= 0);

  // double s0bar = (s + s0)/2;
  // double ds = s1 - s0;
  // double alpha_sq = std::pow((u0 - u1)/h, 2);

  // TODO: handle case where sbar0^2 == alpha^2 or prove that it cannot happen

  // double a = s0bar*s0bar - alpha_sq;
  // double b = s0bar*ds;
  // double c = ds/4 - alpha_sq;
  // double disc = b*b - 4*a*c;
  // assert(disc >= 0);
  // double tmp1 = -b/(2*a), tmp2 = std::sqrt(disc)/(2*a);

  printf("warning! not implemented!\n");

  return std::numeric_limits<double>::infinity();
}

double mp1_adj(double u0, double u1, double s0, double s1, double h) {
  check_params(u0, u1, h);
  assert(s0 >= 0);
  assert(s1 >= 0);

  double alpha_sq = std::pow((u0 - u1)/h, 2);
  double s0_sq = s0*s0;
  double ds = s1 - s0;
  double ds_sq = ds*ds;
  double a[] = {
    s0_sq - alpha_sq - 2*s0*ds + ds_sq,
    -4*s0_sq + 2*alpha_sq + 10*s0*ds - 6*ds_sq,
    4*s0_sq - 2*alpha_sq - 20*s0*ds + 17*ds_sq,
    16*s0*ds - 24*ds_sq,
    16*ds_sq
  };

  double lam = -1;  double z[8];
  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(5);
  gsl_poly_complex_solve(a, 5, w, z);
  gsl_poly_complex_workspace_free(w);

  for (int i = 0; i < 4; ++i) {
    if (z[2*i + 1] != 0 || z[2*i] < 0 || z[2*i] > 1) {
      continue;
    }
    lam = z[2*i];
    double lhs = (u0 - u1)*std::sqrt(1 - 2*lam + 2*lam*lam)/h;
    double rhs = -s0*(4*lam*lam - 5*lam + 2) + s1*(4*lam*lam - 3*lam + 1);
    if (std::fabs(lhs - rhs)/std::fabs(lhs) < 1e-6) {
      break;
    }
  }
  return lam == -1 ?
    std::numeric_limits<double>::infinity() :
    (1 - lam)*u0+ lam*u1 + ((1 - lam)*s0 + lam*s1)*h*std::sqrt(1 - 2*lam + 2*lam*lam + 1);
}

double mp1_diag(double u0, double u1, double s0, double s1, double h) {
  check_params(u0, u1, h);
  assert(s0 >= 0);
  assert(s1 >= 0);

  double alpha = std::fabs((u0 - u1)/h);
  double ds = s1 - s0;
  double a[] = {
    (ds - alpha)*(ds + alpha),
    2*s0*ds,
    4*ds*ds + (s0 - alpha)*(s0 + alpha),
    4*s0*ds,
    4*ds*ds
  };

  double z[8];
  gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(5);
  gsl_poly_complex_solve(a, 5, w, z);
  gsl_poly_complex_workspace_free(w);

  double lam = -1;
  for (int i = 0; i < 4; ++i) {
    if (z[2*i + 1] != 0 || z[2*i] < 0 || z[2*i] > 1) {
      continue;
    }
    lam = z[2*i];
    double lhs = (u0 - u1)*std::sqrt(lam*lam + 1)/h;
    double rhs = s1 + 2*s1*lam*lam - s0*(1 - lam + 2*lam*lam);
    if (std::fabs(lhs - rhs)/std::fabs(lhs) < 1e-6) {
      break;
    }
  }
  return lam == -1 ?
    std::numeric_limits<double>::infinity() :
    (1 - lam)*u0+ lam*u1 + ((1 - lam)*s0 + lam*s1)*h*std::sqrt(lam*lam + 1);
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
                double left, double right) {
  int nroots = sturm(polys, left, right);
  if (nroots == 1) {
    double h = 0.1, x0;
    bool foundroot = false;
    x0 = secant(polys[0], left, left + (right - left)*h, left, right, foundroot);
    if (!foundroot) {
      x0 = secant(polys[0], right, right + (left - right)*h, left, right, foundroot);
    }
    assert(foundroot);
    roots[root++] = x0;
  } else if (nroots > 1) {
    double mid = (left + right)/2;
    rec(polys, roots, root, left, mid);
    rec(polys, roots, root, mid, right);
  }
}

void find_quartic_roots(double * a, double * roots, double left, double right) {
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
  rec(polys, roots, root, left, right);

  // printf("roots = {");
  // if (root > 0) {
  //   for (int i = 0; i < root - 1; ++i) {
  //     printf("%g, ", roots[i]);
  //   }
  //   printf("%g}\n", roots[root - 1]);
  // } else {
  //   printf("}\n");
  // }

  for (; root < 5; ++root) {
    roots[root] = -1;
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
