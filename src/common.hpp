#pragma once

// TODO: try to remove these
#include <limits>
#include <type_traits>

#include <assert.h>

template <bool value>
using bool_t = std::integral_constant<bool, value>;

constexpr double pi = 3.141592653589793;

constexpr double int_sqrt(int i) {
  assert(i < 7);
  double _lut[] = {
    0.0,                // i = 0
    1.0,                // i = 1
    1.4142135623730951, // i = 2
    1.7320508075688772, // i = 3
    2.0,                // i = 4
    2.2360679774997897, // i = 5
    2.4494897427831779  // i = 6
  };
  return _lut[i];
}

template <class T>
constexpr T eps = std::numeric_limits<T>::epsilon();

template <class T>
constexpr T inf = std::numeric_limits<T>::infinity();

#ifdef OLIM_DEBUG
#  define OLIM_PROTECTED public
#  define OLIM_PRIVATE public
#else
#  define OLIM_PROTECTED protected
#  define OLIM_PRIVATE private
#endif
