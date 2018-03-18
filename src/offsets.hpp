#ifndef __OFFSETS_HPP__
#define __OFFSETS_HPP__

#pragma GCC diagnostic ignored "-Wunused-variable"

namespace {
  template <int d>
  struct offsets_size {
    static_assert(d == 2 || d == 3, "invalid number of dimensions");
    static int size;
  };
}

namespace {
  template <int d> constexpr int di[offsets_size<d>::size];
  template <int d> constexpr int dj[offsets_size<d>::size];
  template <int d> constexpr int dk[offsets_size<d>::size];
}

#include "offsets.impl.hpp"

#pragma GCC diagnostic pop

#endif // __OFFSETS_HPP__
