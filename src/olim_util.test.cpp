#include "olim_util.hpp"
#include "test.hpp"

void rhr_3d_22_works() {
  using std::sqrt;
  double h = 0.1;
  IS_APPROX_EQUAL(rhr_3d_22(0.1, 0.0, 1.0, h), h*sqrt(2));
  IS_APPROX_EQUAL(rhr_3d_22(0.0, 0.1, 1.0, h), h*sqrt(2));
  IS_APPROX_EQUAL(rhr_3d_22(0.0, 0.0, 1.0, h), h*sqrt(6)/2);
}

int main() {
  rhr_3d_22_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
