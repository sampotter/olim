#pragma once

#include "../vec.hpp"

#ifndef __clang__
#    pragma GCC diagnostic ignored "-Wunused-variable"
#endif

namespace {

template <int n> constexpr int offset_size() {
  if (n == 2) {
    return 8;
  } else if (n == 3) {
    return 26;
  }
}

/**
 * This array stores offset indices into neighborhoods, where each the
 * max-norm cardinality of each offset is 1. The first index is the
 * dimension of the ambient space, and the offsets are partially
 * sorted by increasing nonzero Hamming weight, and there are 3^n - 1
 * offsets total. The ordering by Hamming weight sorts the offsets to
 * induce symmetric neighborhoods of different size; i.e., for n = 3,
 * there are 6 Hamming weight 1 offsets and 12 Hamming weight 2
 * offsets. We can restrict ourselves to the first 18 offsets to use a
 * symmetric neighborhood of size 18 (where the complete max-norm unit
 * ball has 26 offsets).
 */
template <int n> constexpr int offsets[n][offset_size<n>];

/**
 * The first 4 offsets have Hamming weight 1 and form the olim4
 * neighborhood. The next 4 offsets have Hamming weight 2. All offsets
 * together form the olim8 neighborhood.
 */
template <> int offsets<2>[2][offset_size<2>()] = {
  {-1, 0, 1,  0, -1, 1,  1, -1},
  { 0, 1, 0, -1,  1, 1, -1, -1}
};

/**
 * The first 6 offsets have Hamming weight 1 (the 6 "cardinal"
 * neighbors). These correspond to the olim6 neighbors. The next 12
 * have Hamming weight 2. These combined with the Hamming weight 1
 * offsets comprise the olim18 neighborhood. The last 8 have Hamming
 * weight 3. All 26 offsets form the olim26 neighborhood.
 */
template <> int offsets<3>[3][offset_size<3>()] = {
  {1,  0,  0, -1,  0,  0,
   1,  0, -1,  0,  1, -1, -1,  1,  1,  0, -1,  0,
   1, -1, -1,  1,  1, -1, -1,  1},
  {0,  1,  0,  0, -1,  0,
   0,  1,  0, -1,  1,  1, -1, -1,  0,  1,  0, -1,
   1,  1, -1, -1,  1,  1, -1, -1},
  {0,  0,  1,  0,  0, -1,
   1,  1,  1,  1,  0,  0,  0,  0, -1, -1, -1, -1,
   1,  1,  1,  1, -1, -1, -1, -1}
};

template <int n>
inline vec<double, n> get_p(int i) {
  vec<double, n> p;
  for (int j = 0; j < n; ++j) {
    p[j] = offsets<n>[j][i];
  }
  return p;
}

template <int n>
inline vec<int, n> get_offset(int i) {
  vec<int, n> offset;
  for (int j = 0; j < n; ++j) {
    offset[j] = offsets<n>[j][i];
  }
  return offset;
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
