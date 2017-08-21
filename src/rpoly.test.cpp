#include "rpoly.hpp"
#include "test.hpp"

int main() {
  double gtzeros[] = {0.2, 0.4, 0.6, 0.8};
  double a[] = {0.0384, -0.4, 1.4, -2, 1};
  std::complex<double> roots[4];
  int degree = 4;
  test::is_true(!rpoly(a, degree, roots));
  for (int i = 0; i < degree; ++i) {
    test::is_approx_equal(roots[i].real(), gtzeros[i], 8e-15);
    test::is_approx_equal(roots[i].imag(), 0.0, 1e-15);
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
