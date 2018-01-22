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
      T = min(T, TRI(i, j, p0, p1));            \
    }                                           \
  } while (0)

template <class node, class line_updates, class tri_updates, bool adj_updates,
          bool diag_updates>
void olim<node, line_updates, tri_updates, adj_updates,
          diag_updates>::update_impl(int i, int j, double & T)
{
  using std::min;

#ifdef PRINT_UPDATES
  printf("olim::update_impl(i = %d, j = %d)\n", i, j);
#endif

  abstract_node * nb[num_neighbors];
  memset(nb, 0x0, num_neighbors*sizeof(abstract_node *));
  this->get_valid_neighbors(i, j, nb);

  double h = this->get_h(), s = this->speed(i, j), s_[num_neighbors];
  for (int k = 0; k < num_neighbors; ++k) {
    if (nb[k]) {
      s_[k] = this->speed(i + __di(k), j + __dj(k));
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

#ifdef PRINT_UPDATES
  printf("olim::update_impl: T <- %g\n", T);
#endif
}

template <class node, class line_updates, class tri_updates, bool adj_updates,
		  bool diag_updates>
void olim<node, line_updates, tri_updates, adj_updates, diag_updates>::update_impl(
  int i, int j, int src, double & T)
{
  update_impl(i, j, src, T, eikonal::bool_t<diag_updates> {});
}

template <class node, class line_updates, class tri_updates, bool adj_updates,
		  bool diag_updates>
void olim<node, line_updates, tri_updates, adj_updates, diag_updates>::update_impl(
  int i, int j, int src, double & T, eikonal::bool_t<false>)
{
  using std::min;
#ifdef PRINT_UPDATES
  printf("olim::update_impl(i = %d, j = %d, src = %d)\n", i, j, src);
#endif

  int k[3] = {
    (src + num_neighbors - 1) % num_neighbors,
    src,
    (src + 1) % num_neighbors
  };

  abstract_node * nb[3];
  for (int l = 0, a, b; l < 3; ++l) {
    a = i + __di(k[l]), b = j + __dj(k[l]);
    nb[l] = this->in_bounds(a, b) && this->is_valid(a, b) ?
      &this->operator()(a, b) : nullptr;
  }
  assert(nb[1]); // src node is always in bounds and valid

  double h = this->get_h(), s = this->speed(i, j), s_[3];
  for (int l = 0; l < 3; ++l) {
    if (nb[l]) {
      s_[l] = this->speed(i + __di(k[l]), j + __dj(k[l]));
    }
  }

  DO_LINE(1, 1);
  if (nb[0]) T = min(T, min(LINE(0, 1), TRI(0, 1, 01, 10)));
  if (nb[2]) T = min(T, min(LINE(2, 1), TRI(1, 2, 10, 01)));

#ifdef PRINT_UPDATES
  printf("olim::update_impl: T <- %g\n", T);
#endif
}

template <class node, class line_updates, class tri_updates, bool adj_updates,
		  bool diag_updates>
void olim<node, line_updates, tri_updates, adj_updates, diag_updates>::update_impl(
  int i, int j, int src, double & T, eikonal::bool_t<true>)
{
  using std::min;

#ifdef PRINT_UPDATES
  printf("olim::update_impl(i = %d, j = %d, src = %d)\n", i, j, src);
#endif

  int k[5];
  k[2] = src;
  if (src < 4) {
    k[0] = (src + 3) % 4;
    k[1] = src == 0 ? 7 : src + 3;
    k[3] = src + 4;
    k[4] = (src + 1) % 4;
  } else {
    k[0] = ((src - 1) % 4) + 4;
    k[1] = src - 4;
    k[3] = src - 3;
    k[4] = ((src - 3) % 4) + 4;
  }

  abstract_node * nb[5];
  for (int l = 0, a, b; l < 5; ++l) {
    a = i + __di(k[l]), b = j + __dj(k[l]);
    nb[l] = this->in_bounds(a, b) && this->is_valid(a, b) ?
      &this->operator()(a, b) : nullptr;
  }
  assert(nb[2]); // src node is always in bounds and valid

  double h = this->get_h(), s = this->speed(i, j), s_[5];
  for (int l = 0; l < 5; ++l) {
    if (nb[l]) {
      s_[l] = this->speed(i + __di(k[l]), j + __dj(k[l]));
    }
  }

  if (src < 4) {
    DO_LINE(0, 1);
    DO_LINE(1, 2);
    DO_LINE(2, 1);
    DO_LINE(3, 2);
    DO_LINE(4, 1);
    DO_TRI(0, 1, 01, 11);
    DO_TRI(1, 2, 11, 10);
    DO_TRI(2, 3, 01, 11);
    DO_TRI(3, 4, 11, 10);
  } else {
    DO_LINE(0, 2);
    DO_LINE(1, 1);
    DO_LINE(2, 2);
    DO_LINE(3, 1);
    DO_LINE(4, 2);
    DO_TRI(0, 1, 11, 10);
    DO_TRI(1, 2, 01, 11);
    DO_TRI(2, 3, 11, 10);
    DO_TRI(3, 4, 01, 11);
  }

#ifdef PRINT_UPDATES
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
