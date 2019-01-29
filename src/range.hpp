#pragma once

#include "common.hpp"
#include "vec.hpp"

template <int n>
struct range
{
  range(vec<int, n> dims):
    _dims {dims},
    _end_lin {_dims.product()},
    _end {to_vector_index(_end_lin, _dims)}
  {}

  inline range<n> & operator++() {
    ++_lin;
    return *this;
  }

  inline range<n> operator++(int) {
    range<n> tmp {*this};
    ++(*this);
    return tmp;
  }

  inline bool operator!=(range<n> const & other) const {
    return _lin != other._lin || _dims != other._dims;
  }

  inline vec<int, n> operator*() const {
    return to_vector_index(_lin, _dims);
  }

  inline range<n> begin() const {
    return {_dims};
  }

  inline range<n> end() const {
    range<n> tmp {_dims};
    tmp._lin = tmp._end_lin;
    return tmp;
  }
  
OLIM_PRIVATE:
  int _lin {0};
  vec<int, n> _dims;
  int _end_lin;
  vec<int, n> _end;
};
