#ifndef __COMMON_MACROS_HPP__
#define __COMMON_MACROS_HPP__

#include <limits>

#define EPS(T) std::numeric_limits<T>::epsilon()
#define INF(T) std::numeric_limits<T>::infinity()
#define VAL(i) (nb[i]->get_value())

#define GET_MACRO_NAME_3(_1, _2, _3, NAME, ...) NAME

#define SPEED_ARGS(...)                         \
  GET_MACRO_NAME_3(                             \
    __VA_ARGS__,                                \
    SPEED_ARGS_3,                               \
    SPEED_ARGS_2,                               \
    SPEED_ARGS_1)(__VA_ARGS__)

#define SPEED_ARGS_1(i) s, s_[i]
#define SPEED_ARGS_2(i, j) s, s_[i], s_[j]
#define SPEED_ARGS_3(i, j, k) s, s_[i], s_[j], s_[k]

#ifdef EIKONAL_DEBUG
#  define EIKONAL_PROTECTED public
#  define EIKONAL_PRIVATE public
#else
#  define EIKONAL_PROTECTED protected
#  define EIKONAL_PRIVATE private
#endif

#endif // __COMMON_MACROS_HPP__
