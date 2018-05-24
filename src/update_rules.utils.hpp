#ifndef __UPDATE_RULES_UTILS_HPP__
#define __UPDATE_RULES_UTILS_HPP__

#include <cassert>
#include <utility>

template <int dim>
double norm2(double const * p);

template <>
inline double norm2<2>(double const * p) {
  return std::sqrt(p[0]*p[0] + p[1]*p[1]);
}

template <>
inline double norm2<3>(double const * p) {
  return std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

/**
 * This typedef allows us to optionally return whether an update was
 * degenerate or not (which only makes sense for d > 1 updates). This
 * is only done if the COLLECT_STATS compiler flag is set in the CMake
 * build.
 */
template <int d>
struct update_info {
  double value {INF(double)};
  double lambda[d];
#ifdef COLLECT_STATS
  inline bool is_degenerate() const {
    bool deg = false;
    double lam0 = 1;
    for (int i = 0; i < d; ++i) {
      deg = deg || lambda[i] < EPS(double);
      lam0 -= lambda[i];
    }
    return deg || lam0 > 1 - EPS(double);
  }
  bool hierarchical {false};
#endif
};

template <char c>
using ffvec = std::integral_constant<char, c>;

/**
 * Number of bits that are set for integers 0 through 7.
 */
constexpr char nbits[8] = {0, 1, 1, 2, 1, 2, 2, 3};

/**
 * Gets the ith component of the bit vector p.
 */
constexpr char component(char p, char i) {
  // TODO: handle signedness of p
  return (p & 1 << i) >> i;
}

/**
 * Computes Hamming norm of p, interpreted as a bit vector. Assumes
 * that 0 <= p < 8. Implicitly assumes the dimension of the vector is
 * between 1 and 3. For vectors with dimension 1 or 2, make sure the
 * more significant bits are set to zero (i.e., valid vectors of
 * dimension 1 and 2 are 0, 1, 2, and 3).
 */
constexpr char weight(char p) {
  assert(0 <= p && p < 8);
  return nbits[static_cast<int>(p)];
}

/**
 * Compute dot product between p and q, interpreting them as bit
 * vectors. Assumes that 0 <= p < 8 and 0 <= q < 8.
 */
constexpr char dot(char p, char q) {
  assert(0 <= p && p < 8);
  assert(0 <= q && q < 8);
  return nbits[p & q];
}

/**
 * Project coordinates (lam1, lam2) onto nearest barycentric
 * coordinates; i.e., orthogonally project (lam1, lam2) onto the
 * simplex {(x, y): x >= 0, y >= 0, x + y <= 1}.
 */
inline void proj_bary(double & lam1, double & lam2) {
  double const sum = lam1 + lam2;
  if (sum > 1) {
    double const shift = (1 - sum)/2;
    lam1 += shift;
    lam2 += shift;
  }
  lam1 = lam1 < 0 ? 0 : (lam1 > 1 ? 1 : lam1);
  lam2 = lam2 < 0 ? 0 : (lam2 > 1 ? 1 : lam2);
}

#endif // __UPDATE_RULES_UTILS_HPP__
