#ifndef __TEST_HPP__
#define __TEST_HPP__

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <type_traits>

namespace test {
  void is_true(bool prop) {
    if (!prop) {
      fprintf(stdout, "failure: is_true\n");
    }
  }

  template <class T>
  void is_approx_equal(T t, T t_hat, T tol = 1e-7) {
    static_assert(std::is_floating_point<T>::value);
    T const val = std::fabs(t - t_hat)/std::fabs(t);
    if (val > tol) {
      fprintf(stdout, "failure: %g > %g\n", val, tol);
    }
  }
}

#endif // __TEST_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
