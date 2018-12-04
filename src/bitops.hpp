#ifndef __BITOPS_HPP__
#define __BITOPS_HPP__

#include <tuple>

namespace bitops {

template <int n>
using dim = std::integral_constant<int, n>;

template <int i, int... ps>
constexpr int get() {
  return std::get<i>(std::forward_as_tuple(ps...));
}

template <int i, int p>
constexpr int bit() {
  return (p & (1 << i)) >> i;
}

template <int p, int q>
constexpr int p_dot_q(dim<3>) {
  return bit<0, p>()*bit<0, q>() + bit<1, p>()*bit<1, q>() +
    bit<2, p>()*bit<2, q>();
}

template <int p0, int p1, int p2, int i>
constexpr int dPt_dP(dim<3> n) {
  static_assert(0 <= i && i < 3, "need 0 <= i < d");
  if (i == 0) {
    return p_dot_q<p1, p1>(n) - 2*p_dot_q<p1, p0>(n) + p_dot_q<p0, p0>(n);
  } else if (i == 1) {
    return p_dot_q<p1, p2>(n) - p_dot_q<p1, p0>(n) - p_dot_q<p2, p0>(n) + p_dot_q<p0, p0>(n);
  } else {
    return p_dot_q<p2, p2>(n) - 2*p_dot_q<p2, p0>(n) + p_dot_q<p0, p0>(n);
  }
}

}

#endif // __BITOPS_HPP__
