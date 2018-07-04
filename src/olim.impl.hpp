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

#define LINE_FAC(k)                                                     \
  this->template line<2>(VAL(k), SPEED_ARGS(k), h, p0, p_fac, s_fac)

#define DO_LINE_FAC(k) do {                     \
    if (nb[k]) {                                \
      p0[0] = __di(k);                          \
      p0[1] = __dj(k);                          \
      T = min(T, LINE_FAC(k));                  \
    }                                           \
  } while (0)

#define TRI_FAC(k, l)                           \
  this->template tri<2>(                        \
    VAL(k),                                     \
    VAL(l),                                     \
    SPEED_ARGS(k, l),                           \
    h,                                          \
    p0,                                         \
    p1,                                         \
    p_fac,                                      \
    s_fac)                                      \

#define DO_TRI_FAC(a, b) do {                   \
    if (nb[a] && nb[b]) {                       \
      p0[0] = __di(a);                          \
      p0[1] = __dj(a);                          \
      p1[0] = __di(b);                          \
      p1[1] = __dj(b);                          \
      auto tmp = TRI_FAC(a, b);                 \
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
update_impl(node * n, node ** nb, double & T)
{
  using std::min;

  int i_hat = n->get_i(), j_hat = n->get_j();
#if PRINT_UPDATES
  printf("olim::update_impl(i = %d, j = %d)\n", i, j);
#endif

  double h = this->get_h(), s = this->get_speed(i_hat, j_hat), s_[nneib];
  for (int k = 0; k < nneib; ++k) {
    if (nb[k]) {
      s_[k] = this->get_speed(i_hat + __di(k), j_hat + __dj(k));
    }
  }

  if (n->has_parent()) {
    auto n_fac = static_cast<node *>(n->get_parent());
    int i_fac = n_fac->get_i(), j_fac = n_fac->get_j();
    double s_fac = this->get_speed(i_fac, j_fac);
    double p0[2], p1[2];
    double p_fac[2] = {(double) (i_fac - i_hat), (double) (j_fac - j_hat)};

    for (int a = 0, b = 1; a < 4; b = (++a + 1) % 4) {
      DO_LINE(a, 1);
      DO_TRI_FAC(a, b);
    }
    if (diag_updates) {
      for (int a = 4, b = 0, c = 1; a < 8; ++a, c = (++b + 1) % 4) {
        DO_LINE(a, 2);
        DO_TRI_FAC(a, b);
        DO_TRI_FAC(a, c);
      }
    }
  }
  else {
    for (int a = 0, b = 1; a < 4; b = (++a + 1) % 4) {
      DO_LINE(a, 1);
      DO_TRI(a, b, 01, 10);
    }
    if (diag_updates) {
      for (int a = 4, b = 0, c = 1; a < 8; ++a, c = (++b + 1) % 4) {
        DO_LINE(a, 2);
        DO_TRI(a, b, 11, 01);
        DO_TRI(a, c, 11, 10);
      }
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
