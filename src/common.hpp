#pragma once

// TODO: try to remove these
#include <limits>
#include <type_traits>

template <bool value>
using bool_t = std::integral_constant<bool, value>;

constexpr double pi = 3.141592653589793;
constexpr double sqrt2 = 1.4142135623730951;
constexpr double sqrt3 = 1.7320508075688772;
constexpr double sqrt6 = 2.4494897427831779;

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
