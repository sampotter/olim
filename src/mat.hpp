#pragma once

#include "vec.hpp"

// TODO: all this stuff is going to be slow... but hey

// for now, this matrix is going to be stored with row major ordering

template <class T, int m, int n>
struct mat
{
  mat() {}

  mat(mat<T, m, n> const & A) {
    for (int i = 0; i < m*n; ++i) {
      _data[i] = A._data[i];
    }
  }

  mat(T const * t_ptr) {
    for (int i = 0; i < m*n; ++i) {
      _data[i] = t_ptr[i];
    }
  }

  template <class... Ts>
  mat(T const & t, Ts const &... ts) {
    static_assert(sizeof...(ts) == m*n - 1);
    init(t, ts...);
  }

  template <class... Ts>
  inline void init(T const & t, Ts const &... ts) {
    _data[m*n - 1 - sizeof...(ts)] = t;
    init(ts...);
  }

  inline void init() {}

  inline static constexpr mat<T, m, n> zero() {
    mat<T, m, n> A;
    for (int i = 0; i < m*n; ++i) {
      A[i] = 0;
    }
    return A;
  }

  inline T & operator[](int i) {
    return _data[i];
  }

  inline T const & operator[](int i) const {
    return _data[i];
  }

  inline T & operator()(int i, int j) {
    return _data[n*i + j];
  }

  inline T const & operator()(int i, int j) const {
    return _data[n*i + j];
  }

  friend vec<T, m> operator*(mat<T, m, n> const & A, vec<T, n> const & x) {
    auto y = vec<T, m>::zero();
    int offset;
    for (int i = 0; i < m; ++i) {
      offset = n*i;
      for (int j = 0; j < n; ++j) {
        y[i] += A[offset + j]*x[j];
      }
    }
    return y;
  }

  template <int p>
  friend
  mat<T, m, p> operator*(mat<T, m, n> const & A, mat<T, n, p> const & B) {
    auto C = mat<T, m, p>::zero();
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < n; ++k) {
          C(i, j) += A(i, k)*B(k, j);
        }
      }
    }
    return C;
  }

  // TODO: it would be nice to do something like this...
  // vec_view<T, m> column(int j) {
  //   return {&_data[m*j]};
  // }

  T _data[m*n];
};
