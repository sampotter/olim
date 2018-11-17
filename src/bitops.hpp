#ifndef __BITOPS_HPP__
#define __BITOPS_HPP__

#include <tuple>

namespace bitops {

template <int n>
using dim = std::integral_constant<int, n>;

template <int i, char... ps>
constexpr char get() {
  return std::get<i>(std::forward_as_tuple(ps...));
}

template <int i, char p>
constexpr char bit() {
  return (p & (1 << i)) >> i;
}

template <char p, char q>
constexpr char p_dot_q(dim<3>) {
  return bit<0, p>()*bit<0, q>() + bit<1, p>()*bit<1, q>() +
    bit<2, p>()*bit<2, q>();
}

}

#endif // __BITOPS_HPP__
