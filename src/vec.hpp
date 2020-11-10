#pragma once

#include <assert.h>
#include <math.h>

#include "common.hpp"

enum class ordering {
  ROW_MAJOR,
  COLUMN_MAJOR
};

template <class T, int n>
struct vec
{
  vec() {}

  // vec(vec<T, n> const & x) {
  //   for (int i = 0; i < n; ++i) {
  //     _data[i] = x._data[i];
  //   }
  // }



  vec(T const * t_ptr) {
    for (int i = 0; i < n; ++i) {
      _data[i] = t_ptr[i];
    }
  }

  template <class... Ts>
  vec(T const & t, Ts const &... ts) {
    static_assert(sizeof...(ts) == n - 1);
    init(t, ts...);
  }

  template <class... Ts>
  inline void init(T const & t, Ts const &... ts) {
    _data[n - 1 - sizeof...(ts)] = t;
    init(ts...);
  }

  inline void init() {}

  inline static constexpr vec<T, n> one() {
    vec<T, n> x;
    for (int i = 0; i < n; ++i) {
      x[i] = 1;
    }
    return x;
  }

  inline static constexpr vec<T, n> zero() {
    vec<T, n> x;
    for (int i = 0; i < n; ++i) {
      x[i] = 0;
    }
    return x;
  }

  template <class S>
  inline operator vec<S, n>() const {
    vec<S, n> x;
    for (int i = 0; i < n; ++i) {
      x[i] = static_cast<S>(_data[i]);
    }
    return x;
  }

  inline T & operator[](int i) {
    return _data[i];
  }

  inline T const & operator[](int i) const {
    return _data[i];
  }

  inline bool operator==(vec<T, n> const & x) const {
    for (int i = 0; i < n; ++i) {
      if (x[i] != _data[i]) {
        return false;
      }
    }
    return true;
  }

  inline bool operator!=(vec<T, n> const & x) const {
    for (int i = 0; i < n; ++i) {
      if (x[i] != _data[i]) {
        return true;
      }
    }
    return false;
  }

  inline bool operator<(vec<T, n> const & x) const {
    for (int i = 0; i < n; ++i) {
      if (_data[i] >= x[i]) {
        return false;
      }
    }
    return true;
  }

  inline bool operator<=(vec<T, n> const & x) const {
    for (int i = 0; i < n; ++i) {
      if (_data[i] > x[i]) {
        return false;
      }
    }
    return true;
  }

  inline vec<T, n> & operator*=(T const & t) {
    for (int i = 0; i < n; ++i) {
      _data[i] *= t;
    }
    return *this;
  }

  friend vec<T, n> operator*(vec<T, n> x, T const & t) {
    x *= t;
    return x;
  }

  friend vec<T, n> operator*(T const & t, vec<T, n> x) {
    x *= t;
    return x;
  }

  inline vec<T, n> & operator/=(T const & t) {
    for (int i = 0; i < n; ++i) {
      _data[i] /= t;
    }
    return *this;
  }

  friend vec<T, n> operator/(vec<T, n> x, T const & t) {
    x /= t;
    return x;
  }

  friend vec<T, n> operator/(T const & t, vec<T, n> x) {
    x /= t;
    return x;
  }

  friend T operator*(vec<T, n> const & x, vec<T, n> const & y) {
    T t {0};
    for (int i = 0; i < n; ++i) {
      t += x[i]*y[i];
    }
    return t;
  }

  inline vec<T, n> & operator+=(vec<T, n> const & x) {
    for (int i = 0; i < n; ++i) {
      _data[i] += x[i];
    }
    return *this;
  }

  friend vec<T, n> operator+(vec<T, n> x, vec<T, n> const & y) {
    x += y;
    return x;
  }

  inline vec<T, n> & operator-=(vec<T, n> const & x) {
    for (int i = 0; i < n; ++i) {
      _data[i] -= x[i];
    }
    return *this;
  }

  friend vec<T, n> operator-(vec<T, n> x, vec<T, n> const & y) {
    x -= y;
    return x;
  }

  inline T norm1() const {
    T t {0};
    for (int i = 0; i < n; ++i) {
      t += fabs(_data[i]);
    }
    return t;
  }

  inline T norm2sq() const {
    T t {0};
    for (int i = 0; i < n; ++i) {
      t += _data[i]*_data[i];
    }
    return t;
  }

  inline T norm2() const {
    return sqrt(norm2sq());
  }

  inline T normi() const {
    T t {0};
    for (int i = 0; i < n; ++i) {
      t = fmax(t, fabs(_data[i]));
    }
    return t;
  }

  inline vec<T, n> floor() const {
    vec<T, n> x;
    for (int i = 0; i < n; ++i) {
      x[i] = floor(_data[i]);
    }
    return x;
  }

  inline T product() const {
    T t {1};
    for (int i = 0; i < n; ++i) {
      t *= _data[i];
    }
    return t;
  }

OLIM_PRIVATE:
  T _data[n];
};

template <class T>
using vec2 = vec<T, 2>;

template <class T>
using vec3 = vec<T, 3>;

template <class S, class T, int n>
auto operator-(vec<S, n> const & x, vec<T, n> const & y) {
  vec<decltype(x[0] - y[0]), n> z;
  for (int i = 0; i < n; ++i) {
    z[i] = x[i] - y[i];
  }
  return z;
}

template <class T>
vec3<T> operator^(vec3<T> const & x, vec3<T> const & y) {
  return {x[1]*y[2] - x[2]*y[1], x[2]*y[0] - x[0]*y[2], x[0]*y[1] - x[1]*y[0]};
}

template <class T, int n>
T dist1(vec<T, n> const & u, vec<T, n> const & v) {
  return (u - v).norm1();
}

template <class T, int n>
T dist2(vec<T, n> const & u, vec<T, n> const & v) {
  return (u - v).norm2();
}

template <class T, int n>
T dist2sq(vec<T, n> const & u, vec<T, n> const & v) {
  return (u - v).norm2sq();
}

template <class T, int n>
T disti(vec<T, n> const & u, vec<T, n> const & v) {
  return (u - v).normi();
}

// TODO: implement both column-major and row-major ordering correctly

template <ordering order>
inline int to_linear_index(vec<int, 1> inds, vec<int, 1>) {
  return inds[0];
}

template <ordering order>
inline int to_linear_index(vec2<int> inds, vec2<int> dims);

template <>
inline int
to_linear_index<ordering::ROW_MAJOR>(vec2<int> inds, vec2<int> dims)
{
  return inds[1] + dims[1]*inds[0];
}

template <>
inline int
to_linear_index<ordering::COLUMN_MAJOR>(vec2<int> inds, vec2<int> dims)
{
  return inds[0] + dims[0]*inds[1];
}

template <ordering order>
inline int
to_linear_index(vec3<int> inds, vec3<int> dims);

template <>
inline int
to_linear_index<ordering::ROW_MAJOR>(vec3<int> inds, vec3<int> dims)
{
  return inds[2] + dims[2]*(inds[1] + dims[1]*inds[0]);
}

template <>
inline int
to_linear_index<ordering::COLUMN_MAJOR>(vec3<int> inds, vec3<int> dims)
{
  return inds[0] + dims[0]*(inds[1] + dims[1]*inds[2]);
}

template <ordering order>
inline vec<int, 1>
to_vector_index(int lin, vec<int, 1> dims)
{
  return lin % dims[0];
}

template <ordering order>
inline vec2<int> to_vector_index(int lin, vec2<int> dims);

template <>
inline vec2<int>
to_vector_index<ordering::ROW_MAJOR>(int lin, vec2<int> dims)
{
  return {lin/dims[1], lin % dims[1]};
}

template <>
inline vec2<int>
to_vector_index<ordering::COLUMN_MAJOR>(int lin, vec2<int> dims)
{
  return {lin % dims[0], lin/dims[0]};
}

template <ordering order>
inline vec3<int> to_vector_index(int lin, vec3<int> dims);

template <>
inline vec3<int>
to_vector_index<ordering::ROW_MAJOR>(int lin, vec3<int> dims)
{
  return {lin/(dims[2]*dims[1]), lin/dims[2] % dims[1], lin % dims[2]};
}

template <>
inline vec3<int>
to_vector_index<ordering::COLUMN_MAJOR>(int lin, vec3<int> dims)
{
  return {lin % dims[0], lin/dims[0] % dims[1], lin/(dims[0]*dims[1])};
}

// TODO: really need to find a better home for this stuff

/**
 * Number of bits that are set for integers 0 through 7.
 */
constexpr int nbits[8] = {0, 1, 1, 2, 1, 2, 2, 3};

/**
 * Compute dot product between p and q, interpreting them as bit
 * vectors. Assumes that 0 <= p < 8 and 0 <= q < 8.
 */
constexpr int dot(int p, int q) {
  assert(0 <= p && p < 8);
  assert(0 <= q && q < 8);
  return nbits[p & q];
}
