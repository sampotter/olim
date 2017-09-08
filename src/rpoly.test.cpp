#include "rpoly.hpp"
#include "test.hpp"

void rpoly_works_1() {
  double gtzeros[] = {0.2, 0.4, 0.6, 0.8};
  double a[] = {0.0384, -0.4, 1.4, -2, 1};
  std::complex<double> roots[4];
  int degree = 4;
  IS_TRUE(!rpoly(a, degree, roots));
  for (int i = 0; i < degree; ++i) {
    IS_APPROX_EQUAL(roots[i].real(), gtzeros[i], 8e-15);
    IS_APPROX_EQUAL(roots[i].imag(), 0.0, 1e-15);
  }
}

void rpoly_works_2() {
  double a[] = {
    -0.99000016666583334,
    0.0098791278581576019,
    -0.013920963747804107,
    0.019758255716315204,
    0.100099996666711110551
  };
  std::complex<double> roots[4];
  IS_TRUE(!rpoly(a, 4, roots));
}

int main() {
  rpoly_works_1();
  rpoly_works_2();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
