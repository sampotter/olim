#include <gtest/gtest.h>

#include "bitops.hpp"
#include "common.hpp"

#define P001 1
#define P010 2
#define P011 3
#define P100 4
#define P101 5
#define P110 6
#define P111 7

TEST (bitops, Q_works) {
  auto N = bitops::dim<3> {};
  double Q00 = bitops::Q<P001, P010, P100, 0, 0>(N);
  double Q10 = bitops::Q<P001, P010, P100, 1, 0>(N);
  double Q20 = bitops::Q<P001, P010, P100, 2, 0>(N);
  double Q01 = bitops::Q<P001, P010, P100, 0, 1>(N);
  double Q11 = bitops::Q<P001, P010, P100, 1, 1>(N);
  double Q21 = bitops::Q<P001, P010, P100, 2, 1>(N);
  ASSERT_NEAR(Q00, -int_sqrt(2)/2, 2e-16);
  ASSERT_NEAR(Q10, +int_sqrt(2)/2, 2e-16);
  ASSERT_NEAR(Q20, 0, 2e-16);
  ASSERT_NEAR(Q01, -1/int_sqrt(6), 2e-16);
  ASSERT_NEAR(Q11, -1/int_sqrt(6), 2e-16);
  ASSERT_NEAR(Q21, +2/int_sqrt(6), 2e-16);
}

TEST (bitops, R_works) {
  auto N = bitops::dim<3> {};
  double R00 = bitops::R<P001, P010, P100, 0>(N);
  double R01 = bitops::R<P001, P010, P100, 1>(N);
  double R11 = bitops::R<P001, P010, P100, 2>(N);
  ASSERT_NEAR(R00, int_sqrt(2), 2e-16);
  ASSERT_NEAR(R01, int_sqrt(2)/2, 2e-16);
  ASSERT_NEAR(R11, int_sqrt(3)/int_sqrt(2), 2e-16);
}
