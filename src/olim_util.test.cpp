#include "olim_util.hpp"
#include "speed_funcs.hpp"
#include "test.hpp"

void figuring_it_out() {}

void sturm_works() {
  double a[] = {
    -0.0019249999999999998,
    -0.0010000000000000002,
    0.00087500000000000078,
    -0.0020000000000000005,
    0.00040000000000000018
  };
  double b[] = {
    -0.0010000000000000002,
    0.0017500000000000016,
    -0.0060000000000000019,
    0.0016000000000000007
  };
  double c[] = {
    0.0022374999999999999,
    0.00020312499999999986,
    0.0014374999999999998
  };
  double d[] = {
    -0.0086910396975425352,
    -0.0051202079395085099
  };
  double e[] = {
    -0.006034391692523379
  };
  double * polys[] = {a, b, c, d, e};
  test::is_true(sturm(polys, 0, 1) == 0);
}

int main() {
  figuring_it_out();
  sturm_works();
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
