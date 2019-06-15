#include <gtest/gtest.h>

#include "mat.hpp"

TEST (mat, init_works) {
  mat<int, 2, 3> A {
    1, 2, 3,
    4, 5, 6
  };
  ASSERT_EQ(A[0], 1);
  ASSERT_EQ(A[1], 2);
  ASSERT_EQ(A[2], 3);
  ASSERT_EQ(A[3], 4);
  ASSERT_EQ(A[4], 5);
  ASSERT_EQ(A[5], 6);
  ASSERT_EQ(A(0, 0), 1);
  ASSERT_EQ(A(0, 1), 2);
  ASSERT_EQ(A(0, 2), 3);
  ASSERT_EQ(A(1, 0), 4);
  ASSERT_EQ(A(1, 1), 5);
  ASSERT_EQ(A(1, 2), 6);
}

TEST (mat, mat_vec_multiply_works) {
  mat<int, 2, 3> A {
    1, 2, 3,
    4, 5, 6
  };
  vec<int, 3> x {1, 2, 3};
  vec<int, 2> y = A*x;
  ASSERT_EQ(y[0], 14);
  ASSERT_EQ(y[1], 32);
}

TEST (mat, mat_mat_multiply_works) {
  mat<int, 2, 3> A = {
    1, 2, 3,
    4, 5, 6
  };
  mat<int, 3, 2> B = {
    6, 5,
    4, 3,
    2, 1
  };
  mat<int, 2, 2> C = A*B;
  mat<int, 3, 3> D = B*A;
  ASSERT_EQ(C(0, 0), 20);
  ASSERT_EQ(C(1, 0), 56);
  ASSERT_EQ(C(0, 1), 14);
  ASSERT_EQ(C(1, 1), 41);
  ASSERT_EQ(D(0, 0), 26);
  ASSERT_EQ(D(0, 1), 37);
  ASSERT_EQ(D(0, 2), 48);
  ASSERT_EQ(D(1, 0), 16);
  ASSERT_EQ(D(1, 1), 23);
  ASSERT_EQ(D(1, 2), 30);
  ASSERT_EQ(D(2, 0), 6);
  ASSERT_EQ(D(2, 1), 9);
  ASSERT_EQ(D(2, 2), 12);
}
