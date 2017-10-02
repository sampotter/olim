#include "basic_marcher.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "common.macros.hpp"

void basic_marcher::update_impl(int i, int j, double & T) {
  using std::min;

  abstract_node * nb[4] = {nullptr, nullptr, nullptr, nullptr};
  get_valid_neighbors(i, j, nb);
  double sh = get_h()*speed(i, j);

  double T1 = min(nb[0] ? VAL(0) : INF(T1), nb[2] ? VAL(2) : INF(T1));
  double T2 = min(nb[1] ? VAL(1) : INF(T2), nb[3] ? VAL(3) : INF(T2));

  bool T1_inf = std::isinf(T1), T2_inf = std::isinf(T2);

  if (!T1_inf && !T2_inf) {
    double diff = T1 - T2, disc = 2*sh*sh - diff*diff;
    T = disc > 0 ? std::min(T, (T1 + T2 + std::sqrt(disc))/2) : T;
  } else if (T1_inf) {
    T = std::min(T, T2 + sh);
  } else if (T2_inf) {
    T = std::min(T, T1 + sh);
  } else {
    assert(false);
  }
}

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
