#pragma once

#include "../offsets.hpp"
#include "../vec.hpp"

#ifndef __clang__
#    pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace eikonal {

template <int n>
inline vec<double, n> get_p(int i) {
  return detail::get_offset<n>(i);
}

template <int n>
inline int get_linear_offset(vec<int, n> offset);

// TODO: the ordering of the these implementations is mixed up

template <>
inline int get_linear_offset<2>(vec2<int> offset) {
  static constexpr int _lut[] = {7, 0, 4, 3, -1, 1, 6, 2, 5};
  return _lut[3*(offset[0] + 1) + offset[1] + 1];
}

template <>
inline int get_linear_offset<3>(vec3<int> offset) {
  static constexpr int _lut[] = {
    24, 16, 23, 17, 5, 15, 25, 14, 22, 12, 3, 11, 4, -1, 1,
    13, 0, 10, 20, 8, 19, 9, 2, 7, 21, 6, 18
  };
  return _lut[3*3*(offset[2] + 1) + 3*(offset[0] + 1) + offset[1] + 1];
}

}

#ifndef __clang__
#    pragma GCC diagnostic pop
#endif
