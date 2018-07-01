#ifndef __OLIM_IMPL_HPP__
#define __OLIM_IMPL_HPP__

#include <src/config.hpp>

#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"
#include "update_rules.utils.hpp"

#define __di(k) di<2>[k]
#define __dj(k) dj<2>[k]

#define P01 1
#define P10 2
#define P11 3

#define LINE(i, d) this->template line<d>(VAL(i), SPEED_ARGS(i), h)

#define DO_LINE(i, d) do {                      \
    if (nb[i]) {                                \
      T = min(T, LINE(i, d));                   \
    }                                           \
  } while (0)

#define TRI(i, j, p0, p1)                       \
  this->tri(                                    \
    VAL(i),                                     \
    VAL(j),                                     \
    SPEED_ARGS(i, j),                           \
    h,                                          \
    ffvec<P ## p0> {},                          \
    ffvec<P ## p1> {})

#define DO_TRI(i, j, p0, p1) do {               \
    if (nb[i] && nb[j]) {                       \
      auto tmp = TRI(i, j, p0, p1);             \
      T = min(T, tmp.value);                    \
    }                                           \
  } while (0)

template <
  class node,
  class line_updates,
  class tri_updates,
  bool adj_updates,
  bool diag_updates>
void
olim<node, line_updates, tri_updates, adj_updates, diag_updates>::
update_impl(int i, int j, abstract_node ** nb, double & T)
{
  using std::min;

#if PRINT_UPDATES
  printf("olim::update_impl(i = %d, j = %d)\n", i, j);
#endif

  double h = this->get_h(), s = this->get_speed(i, j), s_[nneib];
  for (int k = 0; k < nneib; ++k) {
    if (nb[k]) {
      s_[k] = this->get_speed(i + __di(k), j + __dj(k));
    }
  }

  for (int i = 0, j = 1; i < 4; j = (++i + 1) % 4) {
    DO_LINE(i, 1);
    DO_TRI(i, j, 01, 10);
  }
  if (diag_updates) {
    for (int i = 4, j = 0, k = 1; i < 8; ++i, k = (++j + 1) % 4) {
      DO_LINE(i, 2);
      DO_TRI(i, j, 11, 01);
      DO_TRI(i, k, 11, 10);
    }
  }

#if PRINT_UPDATES
  printf("olim::update_impl: T <- %g\n", T);
#endif
}

#undef __di
#undef __dj

#undef P01
#undef P10
#undef P11
#undef LINE
#undef DO_LINE
#undef TRI

#undef DO_TRI

#endif // __OLIM_IMPL_HPP__
