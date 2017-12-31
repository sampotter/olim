#include "basic_marcher.hpp"

#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"

void basic_marcher::update_impl(int i, int j, double & T) {
  using std::min;

#ifdef PRINT_UPDATES
  printf("basic_marcher::update_impl(i = %d, j = %d)\n", i, j);
#endif

  abstract_node * nb[4] = {nullptr, nullptr, nullptr, nullptr};
  get_valid_neighbors(i, j, nb);
  double sh = get_h()*speed(i, j);

  double T1 = min(nb[0] ? VAL(0) : INF(double), nb[2] ? VAL(2) : INF(double));
  double T2 = min(nb[1] ? VAL(1) : INF(double), nb[3] ? VAL(3) : INF(double));

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

#ifdef PRINT_UPDATES
  printf("basic_marcher::update_impl: T <- %g\n", T);
#endif
}
