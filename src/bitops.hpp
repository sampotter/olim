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

}

#endif // __BITOPS_HPP__
