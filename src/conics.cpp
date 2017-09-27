#include "conics.hpp"

#include <cassert>

#include <armadillo>
#include <gsl/gsl_poly.h>

static arma::mat conic_matrix_from_coefs(double const * Q) {
  arma::mat A(3, 3);
  A(0, 0) = Q[0];
  A(1, 0) = Q[1]/2;
  A(2, 0) = Q[3]/2;
  A(0, 1) = A(1, 0);
  A(1, 1) = Q[2];
  A(2, 1) = Q[4]/2;
  A(0, 2) = A(2, 0);
  A(1, 2) = A(2, 1);
  A(2, 2) = Q[5];
  return A;
}

static arma::mat sym_adjoint(arma::mat const & A) {
  assert(A.n_rows == 3);
  assert(A.n_cols == 3);
  arma::mat B(3, 3);

  double a = A(0, 0), b = A(0, 1), c = A(1, 1), d = A(0, 2),
    e = A(1, 2), f = A(2, 2);
  
  B(0, 0) = c*f - e*e;
  B(0, 1) = -b*f + e*d;
  B(0, 2) = b*e - c*d;
  B(1, 0) = B(0, 1);
  B(1, 1) = a*f - d*d;
  B(1, 2) = -a*e + b*d;
  B(2, 0) = B(0, 2);
  B(2, 1) = B(1, 2);
  B(2, 2) = a*c - b*b;

  return B;
}

static arma::mat cross_matrix(arma::vec const & p) {
  arma::mat A(3, 3, arma::fill::zeros);
  A(0, 1) = p(2);
  A(0, 2) = -p(1);
  A(1, 0) = -p(2);
  A(1, 2) = p(0);
  A(2, 0) = p(1);
  A(2, 1) = -p(0);
  return A;
}

static bool split_deg_conic(arma::mat const & A, arma::vec & m, arma::vec & l) {
  int argmax;
  arma::mat B;

  if (rank(A) == 1) {
    B = A;
  } else {
    B = -sym_adjoint(A);
    argmax = abs(B.diag()).index_max();
    if (B(argmax, argmax) < 0) {
      return false;
    }
    B = A + cross_matrix(B.col(argmax))/sqrt(B(argmax, argmax));
  }

  arma::uvec s = ind2sub(size(B), abs(B).index_max());
  l = B.row(s(0)).t();
  m = B.col(s(1));

  return true;
}

static void get_points_on_line(arma::vec const & l, arma::vec & p1,
                               arma::vec & p2) {
  p1.zeros();
  p2.zeros();
  if (l(0) == 0 && l(1) == 0) {
    p1(0) = 1;
    p2(1) = 1;
  } else {
    p2(0) = -l(1);
    p2(1) = l(0);
    if (fabs(l(0)) < fabs(l(1))) {
      p1(1) = -l(2);
      p1(2) = l(1);
    } else {
      p1(0) = -l(2);
      p1(2) = l(0);
    }
  }
}

static bool intersect_conic_with_line(arma::mat const & A, arma::vec const & l,
                                      arma::vec & p1, arma::vec & p2) {
  get_points_on_line(l, p1, p2);
  double d11 = dot(p1, A*p1), d12 = dot(p1, A*p2), d22 = dot(p2, A*p2);
  double delta = d12*d12 - d11*d22, sqrtdelta = sqrt(delta);
  if (delta < 0) return false;
  double numer = sqrtdelta - d12, denom = d22;
  double t1 = numer == 0 && denom == 0 ? 0 : numer/denom;
  numer = -(sqrtdelta + d12);
  double t2 = numer == 0 && denom == 0 ? 0 : numer/denom;
  arma::vec dp1 = t1*p2, dp2 = t2*p2;
  p2 = p1 + dp1;
  p1 = p1 + dp2;
  return true;
}

bool intersect_conics(double const * Q1, double const * Q2, double * P, int & n) {
  using namespace arma;

  mat const A1 = conic_matrix_from_coefs(Q1);
  mat const A2 = conic_matrix_from_coefs(Q2);
  auto const rank1 = rank(A1), rank2 = rank(A2);

  vec m(3), l(3);

  bool split = false, isect_A1 = true;
  if (rank1 == 3 && rank2 == 3) {
    mat const X = A1*inv(A2);

    double const a = trace(X);
    double const b = X(0, 0)*X(1, 1) - X(0, 1)*X(1, 0) + X(1, 1)*X(2, 2) -
      X(1, 2)*X(2, 1) + X(0, 0)*X(2, 2) - X(0, 2)*X(2, 0);
    double const c = det(X);

    double roots[3];
    int nroots = gsl_poly_solve_cubic(a, b, c, &roots[0], &roots[1], &roots[2]);

    int i = 0;
    do {
      double root = roots[i++];
      split = split_deg_conic(A1 + root*A2, m, l);
    } while (i < nroots && !split);
  } else if (rank1 < rank2) {
    split = split_deg_conic(A1, m, l);
    isect_A1 = false;
  } else {
    split = split_deg_conic(A2, m, l);
  }
  if (!split) {
    return false;
  }

  vec p1(3), p2(3);
  
  n = 0;
  if (intersect_conic_with_line(isect_A1 ? A1 : A2, m, p1, p2)) {
    P[2*n] = p1[0]/p1[2];
    P[2*n++ + 1] = p1[1]/p1[2];
    P[2*n] = p2[0]/p2[2];
    P[2*n++ + 1] = p2[1]/p2[2];
  }
  if (intersect_conic_with_line(isect_A1 ? A1 : A2, l, p1, p2)) {
    P[2*n] = p1[0]/p1[2];
    P[2*n++ + 1] = p1[1]/p1[2];
    P[2*n] = p2[0]/p2[2];
    P[2*n++ + 1] = p2[1]/p2[2];
  }

  return true;
}

bool arma_rootfinder::intersect_conics(double const * Q1, double const * Q2,
                                       double * P, int & n) const {
  return ::intersect_conics(Q1, Q2, P, n);
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
