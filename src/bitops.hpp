#ifndef __BITOPS_HPP__
#define __BITOPS_HPP__

#include <math.h>

// TODO: try to remove this
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

template <int p0, int p1, int p2, int i, int j>
constexpr double Q(dim<3> n)
{
  static_assert(0 <= i && i < 3, "need 0 <= i < d");

  double r[3]={0,0,0}, q1[3]={0,0,0}, q2[3]={0,0,0}, z[3]={0,0,0};

  r[0] = sqrt(p_dot_q<p1, p1>(n) - 2*p_dot_q<p1, p0>(n) + p_dot_q<p0, p0>(n));
  q1[0] = (bit<0, p1>() - bit<0, p0>())/r[0];
  q1[1] = (bit<1, p1>() - bit<1, p0>())/r[0];
  q1[2] = (bit<2, p1>() - bit<2, p0>())/r[0];
  r[1] = q1[0]*(bit<0, p2>() - bit<0, p0>()) +
    q1[1]*(bit<1, p2>() - bit<1, p0>()) + q1[2]*(bit<2, p2>() - bit<2, p0>());
  z[0] = bit<0, p2>() - bit<0, p0>() - r[1]*q1[0];
  z[1] = bit<1, p2>() - bit<1, p0>() - r[1]*q1[1];
  z[2] = bit<2, p2>() - bit<2, p0>() - r[1]*q1[2];
  r[2] = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
  q2[0] = z[0]/r[2];
  q2[1] = z[1]/r[2];
  q2[2] = z[2]/r[2];

  return j == 0 ? q1[i] : q2[i];
}

template <int p0, int p1, int p2, int i>
constexpr double R(dim<3> n)
{
  static_assert(0 <= i && i < 3, "need 0 <= i < d");

  double r[3]={0,0,0}, q1[3]={0,0,0}, q2[3]={0,0,0}, z[3]={0,0,0};

  r[0] = sqrt(p_dot_q<p1, p1>(n) - 2*p_dot_q<p1, p0>(n) + p_dot_q<p0, p0>(n));
  q1[0] = (bit<0, p1>() - bit<0, p0>())/r[0];
  q1[1] = (bit<1, p1>() - bit<1, p0>())/r[0];
  q1[2] = (bit<2, p1>() - bit<2, p0>())/r[0];
  r[1] = q1[0]*(bit<0, p2>() - bit<0, p0>()) +
    q1[1]*(bit<1, p2>() - bit<1, p0>()) + q1[2]*(bit<2, p2>() - bit<2, p0>());
  z[0] = bit<0, p2>() - bit<0, p0>() - r[1]*q1[0];
  z[1] = bit<1, p2>() - bit<1, p0>() - r[1]*q1[1];
  z[2] = bit<2, p2>() - bit<2, p0>() - r[1]*q1[2];
  r[2] = sqrt(z[0]*z[0] + z[1]*z[1] + z[2]*z[2]);
  q2[0] = z[0]/r[2];
  q2[1] = z[1]/r[2];
  q2[2] = z[2]/r[2];

  return r[i];
}

#define __Q(i, j) Q<p0, p1, p2, i, j>(n)

template <int p0, int p1, int p2, int j>
constexpr double Qt_dot_p0(dim<3> n)
{
  int p0_[3] = {bit<0, p0>(), bit<1, p0>(), bit<2, p0>()};
  return __Q(0, j)*p0_[0] + __Q(1, j)*p0_[1] + __Q(2, j)*p0_[2];
}

#undef __Q

#define __Qt_p0(j) Qt_dot_p0<p0, p1, p2, j>(n)

template <int p0, int p1, int p2>
constexpr double exact_soln_numer(dim<3> n)
{
  return p_dot_q<p0, p0>(n) - __Qt_p0(0)*__Qt_p0(0) - __Qt_p0(1)*__Qt_p0(1);
}

#undef __Qt_p0

}

#endif // __BITOPS_HPP__
