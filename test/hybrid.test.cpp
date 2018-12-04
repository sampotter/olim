#include <gtest/gtest.h>

#include <cmath>
#include <tuple>

#include "hybrid.hpp"

double biquad_deriv(double x) {
  return std::pow(x - 0.5, 3) - x + 0.5;
}

TEST (hybrid, biquadratic_test_function) {
  double a = 0.0, b = 1.0, opt;
  hybrid_status status;
  std::tie(opt, status) = hybrid(biquad_deriv, a, b);
  ASSERT_DOUBLE_EQ(opt, 0.5);
  ASSERT_TRUE(status == hybrid_status::OK);
}
