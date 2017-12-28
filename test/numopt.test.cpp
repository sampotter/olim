#include "numopt.hpp"
#include "test.hpp"

#include <cstdio>

using namespace arma;

void test_setdiff() {
  arma::uvec u = {2};
  arma::uvec v = numopt::setdiff(0, 4, u);
  IS_TRUE(v.n_rows == 3);
  IS_TRUE(v(0) == 0);
  IS_TRUE(v(1) == 1);
  IS_TRUE(v(2) == 3);
}

void test_qpez() {
  {
    mat G = {{2, 1}, {1, 2}};
    vec c = {3, 2};
    mat A; A.eye(2, 2);
    vec x = numopt::qpez(G, c, A);
    IS_TRUE(all(abs(x) < 1e-15));
  }
  {
    mat G(2, 2); G.eye(2, 2);
    vec c = {1, 1};
    mat A = {{-1, 1}};
    vec x = numopt::qpez(G, c, A);
    IS_APPROX_EQUAL(x(0), -1.);
    IS_APPROX_EQUAL(x(1), -1.);
  }
  {
    mat G(2, 2); G.eye(2, 2);
    vec c = {-1, -1};
    mat A = {{0, 1}};
    vec x = numopt::qpez(G, c, A);
    IS_APPROX_EQUAL(x(0), 1.);
    IS_APPROX_EQUAL(x(1), 0.);
  }
}

void test_qpi() {
  {
    mat G(2, 2); G.eye(2, 2);
    vec c = {-1./3, -1./3};
    mat A = {{1, 0}, {0, 1}, {-1, -1}};
    vec b = {0, 0, -1};
    vec x = numopt::qpi(G, c, A, b);
    IS_APPROX_EQUAL(x(0), 1./3);
    IS_APPROX_EQUAL(x(1), 1./3);

    c(0) = -1;
    c(1) = -1;
    x = numopt::qpi(G, c, A, b);
    IS_APPROX_EQUAL(x(0), 1./2);
    IS_APPROX_EQUAL(x(1), 1./2);
  }
}

int main() {
  test_setdiff();
  test_qpez();
  test_qpi();
}
