#ifndef __OLIM26_RECT_IMPL_HPP__
#define __OLIM26_RECT_IMPL_HPP__

#include <src/config.hpp>

#include <algorithm>
#if PRINT_UPDATES
#    include <cstdio>
#endif

#include "common.macros.hpp"
#include "olim3d.macros.def.hpp"
#include "olim26.defs.hpp"

#define SPEED_ARGS(...)                         \
  GET_MACRO_NAME_3(                             \
    __VA_ARGS__,                                \
    SPEED_ARGS_3,                               \
    SPEED_ARGS_2,                               \
    SPEED_ARGS_1)(__VA_ARGS__)

#define SPEED_ARGS_1(i) this->estimate_speed(s, s_[i])
#define SPEED_ARGS_2(i, j) this->estimate_speed(s, s_[i], s_[j])
#define SPEED_ARGS_3(i, j, k) this->estimate_speed(s, s_[i], s_[j], s_[k])

// neighbor order:
// degree 1: N, E, U, S, W, D
//           0, 1, 2, 3, 4, 5
// degree 2: UN, UE, US, UW, NE, SE, SW, NW, DN, DE, DS, DW
//           6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17
// degree 3: UNE, USE, USW, UNW, DNE, DSE, DSW, DNW
//           18,  19,  20,  21,  22,  23,  24,  25

template <class node, class update_rules, class speed_estimates>
int olim26_rect<node, update_rules, speed_estimates>::di[] = {
  1, 0, 0, -1, 0, 0,
  1, 0, -1, 0, 1, -1, -1, 1, 1, 0, -1, 0,
  1, -1, -1, 1, 1, -1, -1, 1
};

template <class node, class update_rules, class speed_estimates>
int olim26_rect<node, update_rules, speed_estimates>::dj[] = {
  0, 1, 0, 0, -1, 0,
  0, 1, 0, -1, 1, 1, -1, -1, 0, 1, 0, -1,
  1, 1, -1, -1, 1, 1, -1, -1
};

template <class node, class update_rules, class speed_estimates>
int olim26_rect<node, update_rules, speed_estimates>::dk[] = {
  0, 0, 1, 0, 0, -1,
  1, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, -1,
  1, 1, 1, 1, -1, -1, -1, -1
};

template <class node, class update_rules, class speed_estimates>
void olim26_rect<node, update_rules, speed_estimates>::get_valid_neighbors(
  int i, int j, int k, abstract_node ** nb)
{
  int a, b, c;
  for (int l = 0; l < 26; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->is_valid(a, b, c)) {
      nb[l] = &this->operator()(a, b, c);
    }
  }
}

template <class node, class update_rules, class speed_estimates>
void olim26_rect<node, update_rules, speed_estimates>::stage_neighbors_impl(
  abstract_node * n)
{
  int i = static_cast<node *>(n)->get_i();
  int j = static_cast<node *>(n)->get_j();
  int k = static_cast<node *>(n)->get_k();

#if PRINT_UPDATES
  printf("olim26_rect::stage_neighbors_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  for (int l = 0; l < 26; ++l) {
    this->stage(i + di[l], j + dj[l], k + dk[l]);
  }

  int a, b, c;
  for (int l = 0; l < 26; ++l) {
    a = i + di[l], b = j + dj[l], c = k + dk[l];
    if (this->in_bounds(a, b, c) && !this->operator()(a, b, c).is_valid()) {
      this->update(a, b, c);
    }
  }
}

namespace olim26 {
  int line1tris[6][8] = {
    {UN, UNE, NE, DNE, DN, DNW, NW, UNW}, // N
    {UE, USE, SE, DSE, DE, DNE, NE, UNE}, // E
    {UN, UNE, UE, USE, US, USW, UW, UNW}, // U
    {US, USE, SE, DSE, DS, DSW, SW, USW}, // S
    {UW, USW, SW, DSW, DW, DNW, NW, UNW}, // W
    {DN, DNE, DE, DSE, DS, DSW, DW, DNW}, // D
  };

  int line2tris[12][4] = {
    {U, N, UNE, UNW}, // UN
    {U, E, UNE, USE}, // UE
    {U, S, USE, USW}, // US
    {U, W, UNW, USW}, // UW
    {N, E, UNE, DNE}, // NE
    {S, E, USE, DSE}, // SE
    {S, W, USW, DSW}, // SW
    {N, W, UNW, DNW}, // NW
    {D, N, DNE, DNW}, // DN
    {D, E, DNE, DSE}, // DE
    {D, S, DSE, DSW}, // DS
    {D, W, DNW, DSW}, // DW
  };

  int line3tris[8][6] = {
    {U, UN, N, NE, E, UE}, // UNE
    {U, US, S, SE, E, UE}, // USE
    {U, UW, W, SW, S, US}, // USW
    {U, UN, N, NW, W, UW}, // UNW
    {D, DN, N, NE, E, DE}, // DNE
    {D, DS, S, SE, E, DE}, // DSE
    {D, DS, S, SW, W, DW}, // DSW
    {D, DN, N, NW, W, DW}, // DNW
  };
}

template <class node, class update_rules, class speed_estimates>
void olim26_rect<node, update_rules, speed_estimates>::update_impl(
  int i, int j, int k, double & T)
{
  using namespace olim26;
  using std::min;

#ifdef PRINT_UPDATES
  printf("olim26_rect::update_impl(i = %d, j = %d, k = %d)\n", i, j, k);
#endif

  abstract_node * nb[26];
  memset(nb, 0x0, 26*sizeof(abstract_node *));
  get_valid_neighbors(i, j, k, nb);

  double h = this->get_h(), s = this->speed(i, j, k), s_[26];
  for (int l = 0; l < 26; ++l) {
    if (nb[l]) {
      s_[l] = this->speed(i + di[l], j + dj[l], k + dk[l]);
    }
  }
  
  int l, l0, l1, l2, a, b, * is;

  /**
   * Line updates (degrees 1, 2, and 3)
   */
  for (l = 0; l < 6; ++l) if (nb[l]) LINE1(l);
  for (; l < 18; ++l) if (nb[l]) LINE2(l);
  for (; l < 26; ++l) if (nb[l]) LINE3(l);

  /**
   * Degree (1, 2) and (1, 3) triangle updates
   */
  for (l0 = 0; l0 < 6; ++l0) {
    is = olim26::line1tris[l0];
    if (nb[l0]) {
      for (a = 0, l1 = is[0]; a < 8; l1 = is[++a]) {
        if (nb[l1]) {
          if (a % 2 == 0) TRI12(l0, l1);
          else TRI13(l1, l0);
        }
      }
    }
  }

  /**
   * Degree (2, 3) triangle updates
   */
  for (a = 0, l0 = 18; a < 8; ++a, ++l0) {
    is = olim26::line3tris[a];
    if (nb[l0]) {
      for (b = 3, l1 = is[b]; b < 6; l1 = is[++b]) {
        if (nb[l1]) TRI23(l0, l1);
      }
    }
  }

  /**
   * Degree (1, 2, 3) tetrahedron updates
   */
  for (a = 0, l0 = 18; a < 8; ++a, ++l0) {
    is = olim26::line3tris[a];
    if (nb[l0]) {
      for (b = 0, l1 = is[0], l2 = is[1];
           b < 6;
           ++b, l1 = is[b % 6], l2 = is[(b + 1) % 6]) {
        if (nb[l1] && nb[l2]) {
          if (b % 2 == 0) TETRA123(l1, l2, l0);
          else TETRA123(l2, l1, l0);
        }
      }
    }
  }

#ifdef PRINT_UPDATES
  printf("olim26_rect::update_impl: T <- %g\n", T);
#endif
}

#undef SPEED_ARGS
#undef SPEED_ARGS_1
#undef SPEED_ARGS_2
#undef SPEED_ARGS_3

#include "olim3d.macros.undef.hpp"

#endif // __OLIM26_RECT_IMPL_HPP__

// Local Variables:
// indent-tabs-mode: nil
// c-basic-offset: 2
// End:
