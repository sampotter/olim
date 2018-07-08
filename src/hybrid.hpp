#ifndef __HYBRID_HPP__
#define __HYBRID_HPP__

#include <limits>

enum class hybrid_status {
  OK,
  DEGENERATE
};

template <class F, class T>
std::pair<T, hybrid_status>
hybrid(F const & f, T a, T b, T tol = std::numeric_limits<T>::epsilon());

#include "hybrid.impl.hpp"

#endif // __HYBRID_HPP__
