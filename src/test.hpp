#ifndef __TEST_HPP__
#define __TEST_HPP__

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <type_traits>

#define IS_TRUE(prop) test::is_true(prop, __FILE__, __LINE__)

#define IS_FALSE(prop) test::is_false(prop, __FILE__, __LINE__)

#define IS_APPROX_EQUAL(t, t_hat, ...)                                  \
  test::is_approx_equal(t, t_hat, __FILE__, __LINE__, ##__VA_ARGS__)

namespace test {
  void is_true(bool prop, char const * filename, int line) {
    if (!prop) {
      fprintf(stdout, "failure (%s:%d): is_true\n", filename, line);
      std::abort();
    }
  }

  void is_false(bool prop, char const * filename, int line) {
    if (prop) {
      fprintf(stdout, "failure (%s:%d): is_false\n", filename, line);
      std::abort();
    }
  }

  template <class T>
  void is_approx_equal(T t, T t_hat, char const * filename, int line,
                       T tol = 1e-7) {
    static_assert(std::is_floating_point<T>::value);
    T const val = std::fabs(t - t_hat)/std::fabs(t);
    if (std::isnan(val)) {
      fprintf(stdout, "failure (%s:%d): |%g - %g| == %g > %g\n",
              filename, line, t, t_hat, val, tol);
      std::abort();
    } else if (val > tol) {
      fprintf(stdout, "failure (%s:%d): %g > %g\n", filename, line, val, tol);
      std::abort();
    }
  }
}

#endif // __TEST_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
