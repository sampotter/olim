#ifndef __UPDATE_RULES_TETRA_UPDATES_UTIL_HPP__
#define __UPDATE_RULES_TETRA_UPDATES_UTIL_HPP__

/**
 * Number of bits set for integers 0 through 7.
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
 * that 0 <= p < 8.
 */
constexpr char weight(char p) {
  assert(0 <= p && p < 8);
  return nbits[p];
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

#endif // __UPDATE_RULES_TETRA_UPDATES_UTIL_HPP__
