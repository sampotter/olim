#ifndef __OFFSETS_IMPL_HPP__
#define __OFFSETS_IMPL_HPP__

// neighbor order: N, E, S, W, NE, SE, SW, NW (ordered first by
// degree, then clockwise)

namespace {
  template <> int di<2>[] = {-1, 0, 1, 0, -1, 1, 1, -1};

  template <> int dj<2>[] = {0, 1, 0, -1, 1, 1, -1, -1};

  template <> int di<3>[] = {
    1, 0, 0, -1, 0, 0,
    1, 0, -1, 0, 1, -1, -1, 1, 1, 0, -1, 0,
    1, -1, -1, 1, 1, -1, -1, 1
  };

  template <> int dj<3>[] = {
    0, 1, 0, 0, -1, 0,
    0, 1, 0, -1, 1, 1, -1, -1, 0, 1, 0, -1,
    1, 1, -1, -1, 1, 1, -1, -1
  };

  template <> int dk<3>[] = {
    0, 0, 1, 0, 0, -1,
    1, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, -1,
    1, 1, 1, 1, -1, -1, -1, -1
  };
}

#endif // __OFFSETS_IMPL_HPP__
