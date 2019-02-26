#pragma once

#include "vec.hpp"

template <int n>
struct fac_src
{
  fac_src(vec<double, n> coords, double s):
    coords {coords + decltype(coords)::one()}, s {s} {}

  vec<double, n> coords;
  double s;
};
