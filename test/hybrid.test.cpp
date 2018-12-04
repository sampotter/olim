#include <gtest/gtest.h>

#include <cmath>
#include <tuple>

#include "hybrid.hpp"

/**
 * This test minimizes the function f(x) = (x - 1/2)^4 + (x - 1/2)^2
 * by finding the single zero of f'(x) that lies in [0, 1].
 */
TEST (hybrid, minimize_quartic_test) {
  auto const df = [] (double x) {
    return 4*std::pow(x - 0.5, 3) + 2*(x - 0.5);
  };
  double opt;
  hybrid_status status;
  std::tie(opt, status) = hybrid(df, 0., 1.);
  ASSERT_DOUBLE_EQ(opt, 0.5);
  ASSERT_TRUE(status == hybrid_status::OK);
}

/**
 * This test finds a zero of the quartic f(x) = x0^4 - x0^3 + x0^2 - x0,
 * where x0 = x - 1/3. The only zero of f on [0, 1] is x = 1/3.
 */
TEST (hybrid, find_zero_of_quartic_test) {
  auto const f = [] (double x) {
    using std::pow;
    double x0 = x - 1./3;
    return pow(x0, 4) - pow(x0, 3) + pow(x0, 2) - x0;
  };
  double opt;
  hybrid_status status;
  std::tie(opt, status) = hybrid(f, 0., 1.);
  ASSERT_DOUBLE_EQ(opt, 1./3);
  ASSERT_TRUE(status == hybrid_status::OK);
}
