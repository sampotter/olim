#pragma once

#include <assert.h>

#include "range.hpp"
#include "vec.hpp"

template <class T, int n, ordering ord>
struct grid
{
  using ivec = vec<int, n>;

  grid(ivec dims): _dims {dims}, _data {new T[_dims.product()]} {}

  grid(ivec dims, T const & t): grid {dims} {
    fill(t);
  }

  ~grid() {
    delete[] _data;
  }

  void fill(T const & t) {
    for (int i = 0; i < _dims.product(); ++i) {
      _data[i] = t;
    }
  }

  inline T & operator()(ivec inds) {
    return _data[to_linear_index<ord>(inds, _dims)];
  }

  inline T const & operator()(ivec inds) const {
    return _data[to_linear_index<ord>(inds, _dims)];
  }

  grid<T, n, ord> subgrid(ivec lo, ivec hi) {
    for (int d = 0; d < n; ++d) {
      assert(0 <= lo[d] && lo[d] < hi[d] && hi[d] <= _dims[d]);
    }
    ivec subdims = hi - lo;
    grid g {subdims};
    for (auto ind: range<n, ord> {subdims}) {
      g(ind) = this->operator()(ind + lo);
    }
    return g;
  }

OLIM_PRIVATE:
  ivec _dims;
  T * _data;
};
