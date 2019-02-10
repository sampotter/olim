#pragma once

#include "common.hpp"
#include "vec.hpp"

template <int n, ordering ord = ordering::COLUMN_MAJOR>
struct range
{
  range(vec<int, n> dims):
    _dims {dims},
    _end_lin {_dims.product()},
    _end {to_vector_index<ord>(_end_lin, _dims)}
  {}

  static range<n, ord> cube(int dim) {
    return {dim*vec<int, n>::one()};
  }
  
  inline range<n, ord> & operator++() {
    ++_lin;
    return *this;
  }

  inline range<n, ord> operator++(int) {
    range<n, ord> tmp {*this};
    ++(*this);
    return tmp;
  }

  inline bool operator!=(range<n, ord> const & other) const {
    return _lin != other._lin || _dims != other._dims;
  }

  inline vec<int, n> operator*() const {
    return to_vector_index<ord>(_lin, _dims);
  }

  inline range<n, ord> begin() const {
    return {_dims};
  }

  inline range<n, ord> end() const {
    range<n, ord> tmp {_dims};
    tmp._lin = tmp._end_lin;
    return tmp;
  }
  
OLIM_PRIVATE:
  int _lin {0};
  vec<int, n> _dims;
  int _end_lin;
  vec<int, n> _end;
};
