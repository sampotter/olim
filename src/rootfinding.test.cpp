#include <cmath>

#include "rootfinding.hpp"
#include "test.hpp"

double f(double x) {
  return std::pow(x, 4) - x + 0.3;
}

void test_secant() {
  secant s(f);
  test::is_approx_equal(s.find_root(0.25, 0.26, 1e-15), 0.3091322373030649);
  test::is_approx_equal(s.find_root(0.75, 0.74, 1e-15), 0.8682178749840325);
}

int main() {
  test_secant();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
