#pragma once

#include <array>
#include <math.h>

#include "vec.hpp"

struct no_slow_t {};
constexpr no_slow_t no_slow {};

template <int n>
inline double s0(vec<double, n>) {
  return 1;
}

template <int n>
inline double f0(vec<double, n> x) {
  return x.norm2();
}

template <int n>
inline double s1(vec<double, n> x) {
  return 1 - sin(x.norm2());
}

template <int n>
inline double f1(vec<double, n> x) {
  return cos(x.norm2()) + x.norm2() - 1;
}

template <int n>
inline double s2(vec<double, n> x) {
  return x.norm2();
}

template <int n>
inline double f2(vec<double, n> x) {
  return x.norm2sq()/2;
}

inline double s5(vec2<double> x) {
  return sqrt(pow(x[0], 18) + pow(x[1], 18));
}

inline double f5(vec2<double> x) {
  return (pow(x[0], 10) + pow(x[1], 10))/10;
}

inline double s7(vec2<double> x) {
  x[0] -= 0.900367222589747;
  double aux0 = (x[0] + x[1])/2;
  double aux1 = x[0] + cos(aux0);
  double aux2 = sin(aux0)/2;
  return vec2<double> {{aux1*(1 - aux2), x[1] - aux1*aux2}}.norm2();
}

inline double f7(vec2<double> x) {
  x[0] -= 0.900367222589747;
  return (x[1]*x[1] + pow(x[0] + cos((x[0] + x[1])/2), 2))/2;
}

#ifndef __clang__
#    pragma GCC diagnostic ignored "-Wunused-variable"
#endif

template <int n>
using slow = double(*)(vec<double, n>);

static std::array<slow<2>, 5> slow2s {{s0, s1, s2, s5, s7}};
static std::array<slow<2>, 5> soln2s {{f0, f1, f2, f5, f7}};
static std::array<slow<3>, 3> slow3s {{s0, s1, s2}};
static std::array<slow<3>, 3> soln3s {{f0, f1, s2}};

#ifndef __clang__
#    pragma GCC diagnostic pop
#endif
