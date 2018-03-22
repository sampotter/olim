#ifndef __OFFSETS_HPP__
#define __OFFSETS_HPP__

#ifndef __clang__
#    pragma GCC diagnostic ignored "-Wunused-variable"
#endif

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

#ifndef __clang__
#    pragma GCC diagnostic pop
#endif

#endif // __OFFSETS_HPP__
