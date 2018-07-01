#ifndef __OFFSETS_HPP__
#define __OFFSETS_HPP__

#ifndef __clang__
#    pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace {

template <int d>
constexpr int offset_size();

template <int d> constexpr int di[offset_size<d>()];
template <int d> constexpr int dj[offset_size<d>()];
template <int d> constexpr int dk[offset_size<d>()];

inline int d2l(int di, int dj) {
  static constexpr int _lut[] = {7, 0, 4, 3, -1, 1, 6, 2, 5};
  return _lut[3*(di + 1) + dj + 1];
}

inline int d2l(int di, int dj, int dk) {
  static constexpr int _lut[] = {
    24, 16, 23, 17, 5, 15, 25, 14, 22, 12, 3, 11, 4, -1, 1,
    13, 0, 10, 20, 8, 19, 9, 2, 7, 21, 6, 18
  };
  return _lut[3*3*(dk + 1) + 3*(di + 1) + dj + 1];
}

}

#include "offsets.impl.hpp"

#ifndef __clang__
#    pragma GCC diagnostic pop
#endif

#endif // __OFFSETS_HPP__
