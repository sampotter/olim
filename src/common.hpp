#ifndef __COMMON_DEFS_HPP__
#define __COMMON_DEFS_HPP__

#include <cmath>
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

#ifdef EIKONAL_DEBUG
#  define EIKONAL_PROTECTED public
#  define EIKONAL_PRIVATE public
#else
#  define EIKONAL_PROTECTED protected
#  define EIKONAL_PRIVATE private
#endif

#endif // __COMMON_DEFS_HPP__
