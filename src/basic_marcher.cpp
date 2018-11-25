#include "basic_marcher.hpp"

#include <src/config.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"

void basic_marcher::update_impl(node * n, double & T) {
  using std::min;

#if PRINT_UPDATES
  printf("basic_marcher::update_impl(i = %d, j = %d)\n",n->get_i(),n->get_j());
#else
  (void) n;
#endif

  double sh = this->s_hat*this->get_h();

  double T1 = min(
    nb[0] ? this->nb[0]->get_value() : INF(double),
    nb[2] ? this->nb[2]->get_value() : INF(double));

  double T2 = min(
    nb[1] ? this->nb[1]->get_value() : INF(double),
    nb[3] ? this->nb[3]->get_value() : INF(double));

  bool T1_inf = std::isinf(T1), T2_inf = std::isinf(T2);

  if (!T1_inf && !T2_inf) {
    double diff = T1 - T2, disc = 2*sh*sh - diff*diff;
    T = disc > 0 ? min(T, (T1 + T2 + sqrt(disc))/2) : T;
  } else if (T1_inf) {
    T = min(T, T2 + sh);
  } else if (T2_inf) {
    T = min(T, T1 + sh);
  } else {
    assert(false);
  }

#if PRINT_UPDATES
  printf("basic_marcher::update_impl: T <- %g\n", T);
#endif
}
